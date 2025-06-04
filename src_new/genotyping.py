import functools
import itertools
import math
from typing import Iterator, TypeAlias

import numpy as np
from scipy.stats import binom  # type: ignore


MIN_REP_CNT = 1
MIN_FLANK_LEN = 3
MIN_REP_LEN = 3

Confidences: TypeAlias = tuple[float, float, float, float, float, float, float]


def genotype(
    spanning_observed_counts: list[int], spanning_read_lengths: list[int],
    flanking_observed_counts: list[int], flanking_read_lengths: list[int],
    monoallelic_motif: bool
) -> tuple:
    """This function provides an interface for genotypization step."""
    # 2025-06-04
    # print(f"{spanning_observed_counts=}")
    # print(f"{spanning_read_lengths=}")
    # print(f"{flanking_observed_counts=}")
    # print(f"{flanking_read_lengths=}")
    # print(f"{monoallelic_motif=}")
    # print()

    if len(spanning_observed_counts) == 0 and len(flanking_observed_counts) == 0:
        return (None, ('B', 'B'), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    read_distribution = np.bincount(spanning_read_lengths + flanking_read_lengths, minlength=100)
    flag_spanning_flanking = np.concatenate([
        np.ones_like(spanning_observed_counts, dtype=bool),
        np.zeros_like(flanking_observed_counts, dtype=bool)
    ]).astype(bool)

    model = Inference(read_distribution, None)
    lh_dict = model.real_infer(
        np.array(spanning_observed_counts),
        np.array(flanking_observed_counts),
        monoallelic_motif,
        np.array(spanning_observed_counts + flanking_observed_counts),  # this is stupid
        np.array(spanning_read_lengths + flanking_read_lengths),        # this is stupid as well
        flag_spanning_flanking                                          # ...
    )
    likelihoods, predicted_tmp = model.predict(lh_dict)

    # adjust for no spanning reads (should output Background)
    if len(spanning_observed_counts) == 0:
        predicted_tmp = (0, 0)

    prediction = model.convert_to_sym(predicted_tmp, monoallelic_motif)
    raw_confidence = model.get_confidence(likelihoods, predicted_tmp, monoallelic_motif)
    return likelihoods, prediction, raw_confidence


# All other objects below this line are considered internal a should not be used
class Inference:
    """ Class for inference of alleles. """
    MIN_REPETITIONS = 1
    # default parameters for inference (miSeq default)
    # DEFAULT_MODEL_PARAMS = (0.00716322, 0.000105087, 0.0210812, 0.0001648)
    DEFAULT_MODEL_PARAMS = (0.00716322, 0.000105087, 0.0210812, 0.00716322)  # it actually changes things
    DEFAULT_FIT_FUNCTION = 'linear'

    # TODO: type the members
    def __init__(
        self, read_distribution, params_file,
        str_rep=MIN_REP_CNT, minl_primer1=MIN_FLANK_LEN, minl_primer2=MIN_FLANK_LEN, minl_str=MIN_REP_LEN,
        p_bckg_closed=None, p_bckg_open=None, p_expanded=None
    ):
        """
        Initialization of the Inference class + setup of all models and their probabilities.
        :param read_distribution: ndarray(int) - read distribution
        :param params_file: str - filename of parameters, None for defaults
        :param str_rep: int - length of the STR
        :param minl_primer1: int - minimal length of the left primer
        :param minl_primer2: int - minimal length of the right primer
        :param minl_str: int - minimal length of the STR
        :param p_bckg_closed: float - probability of the background model for closed observation
        :param p_bckg_open: float - probability of the background model for open observation
        :param p_expanded: float - probability of the expanded model (if None it is equal to other models)
        """
        # assign variables
        self.str_rep = str_rep
        self.minl_primer1 = minl_primer1
        self.minl_primer2 = minl_primer2
        self.minl_str = minl_str
        self.read_distribution = read_distribution
        self.sum_reads = np.sum(read_distribution)
        self.params_file = params_file
        self.p_expanded = p_expanded
        self.p_bckg_closed = p_bckg_closed
        self.p_bckg_open = p_bckg_open

    def construct_models(self, min_rep, max_rep, e_model):
        """
        Construct all models needed for current inference.
        :param min_rep: int - minimal allele to model
        :param max_rep: int - maximal allele to model
        :param e_model: int - model for expanded alleles
        :return: None
        """
        # extract params
        model_params = self.DEFAULT_MODEL_PARAMS
        rate_func_str = self.DEFAULT_FIT_FUNCTION
        str_to_func = {'linear': linear_rate, 'const': const_rate, 'exponential': exp_rate, 'square': quadratic_rate}
        rate_func = const_rate
        if rate_func_str in str_to_func.keys():
            rate_func = str_to_func[rate_func_str]

        # save min_rep and max_rep
        self.min_rep = min_rep
        self.max_rep = max_rep  # non-inclusive
        self.max_with_e = e_model + 1  # non-inclusive

        # get models
        mt = model_template(self.max_with_e, model_params, rate_func)
        self.background_model = np.concatenate([
            np.zeros(self.min_rep, dtype=float),
            np.ones(self.max_with_e - self.min_rep, dtype=float) / float(self.max_with_e - self.min_rep)
        ])
        self.expanded_model = mt(self.max_with_e - 1)
        self.allele_models = {i: mt(i) for i in range(min_rep, max_rep)}
        self.models = {'E': self.expanded_model, 'B': self.background_model}
        self.models.update(self.allele_models)

        # get model likelihoods
        open_to_closed = 10.0

        l_others = 1.0
        l_bckg_open = 0.01
        l_exp = 1.01

        l_bckg_model_open = 1.0

        if self.p_expanded is None:
            self.p_expanded = l_exp
        if self.p_bckg_open is None and self.p_bckg_closed is None:
            self.p_bckg_open = l_bckg_open
            self.p_bckg_closed = self.p_bckg_open / open_to_closed
        if self.p_bckg_closed is None:
            self.p_bckg_closed = self.p_bckg_open / open_to_closed
        if self.p_bckg_open is None:
            self.p_bckg_open = self.p_bckg_closed * open_to_closed

        self.model_probabilities = {'E': self.p_expanded, 'B': l_bckg_model_open}
        self.model_probabilities.update({i: l_others for i in self.allele_models.keys()})

    def likelihood_rl(self, rl):
        """
        Likelihood of a read with this length.
        :param rl: int - read length
        :return: float - likelihood of a read this long
        """
        return self.read_distribution[rl] / float(self.sum_reads)

    @staticmethod
    def likelihood_model(model, g):
        """
        Likelihood of a generated allele al from a model of
        :param model: ndarray - model that we evaluate
        :param g: int - observed read count
        :return: float - likelihood of a read coming from this model
        """
        return model[g]

    def likelihood_coverage(self, true_length, rl, closed=True):
        """
        Likelihood of generating a read with this length and this allele.
        :param true_length: int - true number of repetitions of an STR
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return: float - likelihood of a read being generated with this attributes
        """
        whole_inside_str = max(0, true_length * self.str_rep + self.minl_primer1 + self.minl_primer2 - rl + 1)
        # closed_overlapping = max(0, rl - self.minl_primer1 - self.minl_primer2 - true_length * self.str_rep + 1)
        open_overlapping = max(0, rl + true_length * self.str_rep - 2 * self.minl_str + 1)

        assert open_overlapping > whole_inside_str, '%d open %d whole inside %d %d %d' % (
            open_overlapping, whole_inside_str, true_length, rl, self.minl_str)

        return 1.0 / float(open_overlapping - whole_inside_str)

    def likelihood_read_allele(self, model, observed, rl, closed=True):
        """
        Likelihood of generation of read with observed allele count and rl.
        :param model: ndarray - model for the allele
        :param observed: int - observed allele count
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return:
        """
        if closed:
            return (self.likelihood_rl(rl)
                    * self.likelihood_model(model, observed)
                    * self.likelihood_coverage(observed, rl, True))
        number_of_options = 0
        partial_likelihood = 0
        for true_length in itertools.chain(range(observed, self.max_rep), [self.max_with_e - 1]):
            partial_likelihood += (self.likelihood_model(model, true_length)
                                   * self.likelihood_coverage(true_length, rl, False))
            number_of_options += 1

        return self.likelihood_rl(rl) * partial_likelihood / float(number_of_options)

    @functools.lru_cache()
    def likelihood_read(
        self, observed: int, rl: int, model_index1: int, model_index2: int | None = None, closed: bool = True
    ) -> float:
        """
        Compute likelihood of generation of a read from either of those models.
        :param observed: int - observed allele count
        :param rl: int - read length
        :param model_index1: char/int - model index for left allele
        :param model_index2: char/int - model index for right allele or None if mono-allelic
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return: float - likelihood of this read generation
        """
        # TODO: tuto podla mna nemoze byt len tak +, chyba tam korelacia modelov, ale v ramci zjednodusenia asi ok
        allele1_likelihood = (
            self.model_probabilities[model_index1] * self.likelihood_read_allele(self.models[model_index1], observed, rl, closed)
        )
        allele2_likelihood = 0.0 if model_index2 is None else (
            self.model_probabilities[model_index2] * self.likelihood_read_allele(self.models[model_index2], observed, rl, closed)
        )

        p_bckg = self.p_bckg_closed if closed else self.p_bckg_open
        bckgrnd_likelihood = p_bckg * self.likelihood_read_allele(self.models['B'], observed, rl, closed)

        assert not np.isnan(allele2_likelihood)
        assert not np.isnan(allele1_likelihood)
        assert not np.isnan(bckgrnd_likelihood)

        # "tuto" refers to the next line
        return allele1_likelihood + allele2_likelihood + bckgrnd_likelihood

    def real_infer(
        self,
        observed_annots: np.ndarray,
        observed_fa: np.ndarray,
        monoallelic: bool,
        observed_arr: np.ndarray,
        rl_arr: np.ndarray,
        closed_arr: np.ndarray
    ) -> dict[tuple[int | str, int | str], float]:
        """
        Does all the inference,
        computes for which 2 combination of alleles are these annotations and parameters the best.
        argmax_{G1, G2} P(G1, G2 | AL, COV, RL)
            ~ P(AL, COV, RL | G1, G2) * P(G1, G2)
            = prod_{read_i} P(al_i, cov_i, rl_i | G1, G2) * P(G1, G2)
            = independent G1 G2
            = prod_{read_i} P(al_i, cov_i, rl_i | G1) * P(al_i, cov_i, rl_i | G2) * P(G1) * P(G2)
            {here G1, G2 is from possible alleles, background, and expanded, priors are from params}

         P(al_i, cov_i, rl_i | G1) - 2 options:
             1. closed evidence (al_i = X), we know X;
             2. open evidence (al_i >= X), cl_i == True if i is closed

         1.: P(al_i, cov_i, rl_i, cl_i | G1)
            = P(rl_i from read distrib.) * p(allele is al_i | G1) * P(read generated closed evidence | rl_i, al_i)
         2.: P(rl_i is from r.distr.) * P(allele is >= al_i | G1) * P(read generated open evidence | rl_i, al_i)

        :param annotations: list(Annotation) - closed annotated reads (both primers set)
        :param filt_annotations: list(Annotation) - open annotated reads (only one primer set)
        :param index_rep: int - index of a repetition
        :param verbose: bool - print more stuff?
        :param monoallelic: bool - do we have a mono-allelic motif (i.e. chrX/chrY and male sample?)
        :return: dict(tuple(int, int):float) - directory of model indices to their likelihood
        """
        # generate the boundaries:
        overhead = 3
        if len(observed_annots) == 0:
            max_rep = max(observed_fa) + overhead  # non-inclusive
            min_rep = max(self.MIN_REPETITIONS, max(observed_fa) - overhead)  # inclusive
        else:
            max_rep = max(observed_annots) + overhead + 1  # non-inclusive
            min_rep = max(self.MIN_REPETITIONS, min(observed_annots) - overhead)  # inclusive

        # expanded allele
        e_allele = max_rep
        if len(observed_fa) > 0:
            e_allele = max(max_rep, max(observed_fa) + 1)

        # generate all the models
        self.construct_models(min_rep, max_rep, e_allele)

        # go through every model and evaluate:
        evaluated_models = {}
        if monoallelic:
            models = generate_models_one_allele(min_rep, max_rep)
        else:
            models = generate_models(min_rep, max_rep, multiple_backgrounds=True)
        for m1, m2 in models:
            evaluated_models[(m1, m2)] = 0.0
            # go through every read
            for obs, rl, closed in zip(observed_arr, rl_arr, closed_arr):
                lh = self.likelihood_read(obs, rl, m1, None if m2 == 'X' else m2, closed=closed)
                # TODO weighted sum according to the closeness/openness of reads?
                evaluated_models[(m1, m2)] += np.log(lh)

        return evaluated_models

    def predict(self, lh_dict: dict[tuple[int | str, int | str], float]) -> tuple[np.ndarray, tuple[int, int]]:
        # convert to a numpy array:
        lh_array = np.zeros((self.max_rep, self.max_rep + 1))
        for (k1, k2), v in lh_dict.items():
            if k2 == 'X':  # if we have mono-allelic
                k2 = k1
            # B is the smallest, E is the largest!
            if k2 == 'B' or k1 == 'E' or (isinstance(k1, int) and isinstance(k2, int) and k2 < k1):
                k1, k2 = k2, k1
            if k1 == 'B':
                k1 = 0
            if k2 == 'B':
                k2 = 0
            if k1 == 'E':  # only if k2 is 'E' too.
                k1 = 0
            if k2 == 'E':
                k2 = self.max_rep
            lh_array[k1, k2] = v

        # get minimal and maximal likelihood
        ind_good = (lh_array < 0.0) & (lh_array > -1e10) & (lh_array != np.nan)
        if len(lh_array[ind_good]) == 0:
            return lh_array, (0, 0)
        lh_array[~ind_good] = -np.inf
        best = sorted(np.unravel_index(np.argmax(lh_array), lh_array.shape))
        prediction = (int(best[0]), int(best[1]))

        # output best option
        return lh_array, prediction

    def convert_to_sym(self, best: tuple[int, int], monoallelic: bool) -> tuple[int | str, int | str]:
        """
        Convert numeric alleles to their symbolic representations.
        :param best: (int, int) - numeric representation of alleles
        :param monoallelic: bool - if this is monoallelic version
        :return: (int|str, int|str) - symbolic representation of alleles
        """
        # convert it to symbols
        if best[0] == 0 and best[1] == self.max_rep:
            best_sym = ('E', 'E')
        else:
            def fn1(x):
                return 'E' if x == self.max_rep else 'B' if x == 0 else x
            best_sym = tuple(map(fn1, best))

        # if mono-allelic return 'X' as second allele symbol
        if monoallelic:
            best_sym = (best_sym[0], 'X')

        return best_sym

    def get_confidence(
        self, lh_array: np.ndarray, predicted: tuple[int, int], monoallelic: bool = False
    ) -> Confidences:
        """
        Get confidence of a prediction.
        :param lh_array: 2D-ndarray - log likelihoods of the prediction
        :param predicted: tuple(int, int) - predicted alleles
        :param monoallelic: bool - do we have a mono-allelic motif (i.e. chrX/chrY and male sample?)
        :return: tuple[float, float, float | str, float, float, float, float] - prediction confidence of
        all, first, and second allele(s), background and expanded states
        """
        # get confidence
        lh_corr_array = lh_array - np.max(lh_array)
        lh_sum = np.sum(np.exp(lh_corr_array))
        confidence: float = np.exp(lh_corr_array[predicted[0], predicted[1]]) / lh_sum
        confidence1: float
        confidence2: float
        if predicted[0] == predicted[1]:  # same alleles - we compute the probability per allele
            confidence1 = np.sum(np.exp(lh_corr_array[predicted[0], :])) / lh_sum
            confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum
        elif predicted[1] == lh_corr_array.shape[0]:  # expanded allele - expanded is only on one side of the array
            confidence1 = (
                np.sum(np.exp(lh_corr_array[predicted[0], :]))
                + np.sum(np.exp(lh_corr_array[:, predicted[0]]))
                - np.exp(lh_corr_array[predicted[0], predicted[0]])
            ) / lh_sum
            confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum
        else:  # normal behavior - different alleles , no expanded, compute all likelihoods of the alleles
            confidence1 = (
                np.sum(np.exp(lh_corr_array[predicted[0], :]))
                + np.sum(np.exp(lh_corr_array[:, predicted[0]]))
                - np.exp(lh_corr_array[predicted[0], predicted[0]])
            ) / lh_sum
            confidence2 = (
                np.sum(np.exp(lh_corr_array[:, predicted[1]]))
                + np.sum(np.exp(lh_corr_array[predicted[1], :]))
                - np.exp(lh_corr_array[predicted[1], predicted[1]])
            ) / lh_sum

        confidence_back: float = np.exp(lh_corr_array[0, 0]) / lh_sum
        confidence_back_all: float = np.sum(np.exp(lh_corr_array[0, :])) / lh_sum
        confidence_exp: float = np.exp(lh_corr_array[0, self.max_rep]) / lh_sum
        confidence_exp_all: float = np.sum(np.exp(lh_corr_array[:, self.max_rep])) / lh_sum

        if monoallelic:
            # confidence2 = '---'  # TODO: fix this
            confidence2 = np.nan

        return (
            confidence, confidence1, confidence2,
            confidence_back, confidence_back_all,
            confidence_exp, confidence_exp_all
        )


def const_rate(_n, p1=0.0, _p2=1.0, _p3=1.0):
    return p1


def linear_rate(n, p1=0.0, p2=1.0, _p3=1.0):
    return p1 + p2 * n


def quadratic_rate(n, p1=0.0, p2=1.0, p3=1.0):
    return p1 + p2 * n + p3 * n * n


def exp_rate(n, p1=0.0, p2=1.0, p3=1.0):
    return p1 + p2 * math.exp(p3 * n)


def model_template(rng, model_params, rate_func=linear_rate):
    """
    Partial function for model creation.
    :param rng: int - max_range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param rate_func: function - rate function for deletes
    :return: partial function with only 1 parameter - n - target allele number
    """
    return functools.partial(model_full, rng, model_params, rate_func=rate_func)


def generate_models(
    min_rep: int, max_rep: int, multiple_backgrounds: bool = True
) -> Iterator[tuple[int | str, int | str]]:
    """
    Generate all pairs of alleles (models for generation of reads).
    :param min_rep: int - minimal number of repetitions
    :param max_rep: int - maximal number of repetitions
    :param multiple_backgrounds: bool - whether to generate all background states
    :return: generator of allele pairs (numbers or 'E' or 'B')
    """
    for model_index1 in range(min_rep, max_rep):
        for model_index2 in range(model_index1, max_rep):
            yield model_index1, model_index2
        yield model_index1, 'E'
        if multiple_backgrounds:
            yield 'B', model_index1

    yield 'B', 'B'
    yield 'E', 'E'


def generate_models_one_allele(min_rep: int, max_rep: int) -> Iterator[tuple[int | str, int | str]]:
    """
    Generate all pairs of alleles (models for generation of reads).
    :param min_rep: int - minimal number of repetitions
    :param max_rep: int - maximal number of repetitions
    :return: generator of allele pairs (numbers or 'E' or 'B'), 'X' for non-existing allele
    """
    for model_index1 in range(min_rep, max_rep):
        yield model_index1, 'X'

    yield 'B', 'X'
    yield 'E', 'X'


def model_full(rng, model_params, n, rate_func=linear_rate):
    """
    Create binomial model for both deletes and inserts of STRs
    :param rng: int - max_range of distribution
    :param model_params: 4-tuple - parameters for inserts and deletes
    :param n: int - target allele number
    :param rate_func: function - rate function for deletes
    :return: ndarray - combined distribution
    """
    p1, p2, p3, q = model_params
    deletes = binom.pmf(np.arange(rng), n, clip(1 - rate_func(n, p1, p2, p3), 0.0, 1.0))
    inserts = binom.pmf(np.arange(rng), n, q)
    return combine_distribs(deletes, inserts)


def clip(value, minimal, maximal):
    """
    Clips value to range <minimal, maximal>
    :param value: ? - value
    :param minimal: ? - minimal value
    :param maximal: ? - maximal value
    :return: ? - clipped value
    """
    return min(max(minimal, value), maximal)


def combine_distribs(deletes, inserts):
    """
    Combine insert and delete models/distributions
    :param deletes: ndarray - delete distribution
    :param inserts: ndarray - insert distribution
    :return: ndarray - combined array of the same length
    """
    # how much to fill?
    to_fill = sum(deletes == 0.0) + 1
    while to_fill < len(inserts) and inserts[to_fill] > 0.0001:
        to_fill += 1

    # create the end array
    end_distr = np.zeros_like(deletes, dtype=float)

    # fill it!
    for i, a in enumerate(inserts[:to_fill]):
        end_distr[i:] += (deletes * a)[:len(deletes) - i]

    return end_distr
