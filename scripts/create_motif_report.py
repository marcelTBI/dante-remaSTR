from collections import defaultdict
from typing import TypeAlias
from jinja2 import Environment, FileSystemLoader
from typing import Any
import argparse
import json
import textwrap
import os
import sys
import shutil
import pandas as pd


HistReadCounts: TypeAlias = tuple[list[int], list[int], list[int]]
Row: TypeAlias = tuple[str, str | int, str, str | int, str, str, int, int, HistReadCounts]
HeatmapData: TypeAlias = tuple
HistData: TypeAlias = list[int]
PopulationData: TypeAlias = tuple[str, str, int, HeatmapData, HistData, list[Row]]

MIN_GS = 10  # minimal graph size, histograms/heatmaps with less columns look weird


def load_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("")
    )

    arg = parser.add_argument
    arg("-i", "--inputs", nargs="+", help='Path to data.json files (possibly many).')
    arg("-o", "--output_dir", default="./results",
        help='Path to directory where the output will be stored. Default=<input_dir>/result_files')
    arg("-a", "--at-least", type=int, default=0, help="Filter motif sequences with fewer occurences")
    arg("--report_every", default=5, type=int,
        help='Specify how often a progress message should be printed (default=5)')
    arg('-q', '--quiet', action='store_true', help='Don\'t print any progress messages')

    return parser.parse_args()


def collect_jsons(inputs: list[str], output_dir: str) -> dict[str, list[tuple[str, str, dict]]]:
    motif_dict: dict[str, list] = defaultdict(list)
    for filename in inputs:
        print(f"Reading {filename}")
        with open(filename, "r") as f:
            data = json.load(f)

        for motif in data["motifs"]:
            alignment_dir = os.path.dirname(filename) + "/alignments"
            alignment_dir = os.path.relpath(alignment_dir, output_dir)
            motif_dict[motif["motif_id"]].append((data["sample"], alignment_dir, motif))

    return motif_dict


def copy_includes(output_dir: str) -> None:
    include_dir = os.path.dirname(sys.argv[0]) + "/../includes"
    os.makedirs(f'{output_dir}/includes', exist_ok=True)
    shutil.copy2(f'{include_dir}/msa.min.gz.js',            f'{output_dir}/includes/msa.min.gz.js')
    shutil.copy2(f'{include_dir}/plotly-2.14.0.min.js',     f'{output_dir}/includes/plotly-2.14.0.min.js')
    shutil.copy2(f'{include_dir}/jquery-3.6.1.min.js',      f'{output_dir}/includes/jquery-3.6.1.min.js')
    shutil.copy2(f'{include_dir}/datatables.min.js',        f'{output_dir}/includes/datatables.min.js')
    shutil.copy2(f'{include_dir}/styles.css',               f'{output_dir}/includes/styles.css')
    shutil.copy2(f'{include_dir}/w3.css',                   f'{output_dir}/includes/w3.css')
    shutil.copy2(f'{include_dir}/jquery.dataTables.css',    f'{output_dir}/includes/jquery.dataTables.css')


def convert_modules(modules: list) -> list[str]:
    return [f"{seq}[{n}]" for seq, n in modules]


def create_ticks(df: pd.DataFrame) -> list[str]:
    nom1 = set(df.loc[:, ["nomenclature1_len", "nomenclature1_join"]]
                 .itertuples(index=False, name=None))

    nom2 = set(df.loc[:, ["nomenclature2_len", "nomenclature2_join"]]
                 .itertuples(index=False, name=None))

    data = list(nom1.union(nom2))
    data.sort()
    ticktext = [x[1] for x in data]
    return ticktext


def create_main_histrogram(df: pd.DataFrame, ticks: list[str]) -> dict:
    data = [0] * len(ticks)
    counter: set = set()
    for _, row in df.iterrows():
        try:
            idx = ticks.index(row["nomenclature1_join"])
        except ValueError:
            continue
        counter.add(idx)
        data[idx] += 1

    for _, row in df.iterrows():
        try:
            idx = ticks.index(row["nomenclature2_join"])
        except ValueError:
            continue
        counter.add(idx)
        data[idx] += 1

    print(f"Histogram unique nomenclatures: {len(counter)}")
    result = {
        "tickvals": list(range(len(ticks))),
        "ticktext": ticks,
        "x": data
    }
    return result


def create_main_heatmap(df: pd.DataFrame, ticks: list[str]) -> dict:
    data: list[list[int | None]] = [[0 for _ in range(len(ticks))] for _ in range(len(ticks))]
    counter: set = set()
    for _, row in df.iterrows():
        try:
            idx1 = ticks.index(row["nomenclature1_join"])
            idx2 = ticks.index(row["nomenclature2_join"])
        except ValueError:
            continue
        idx1, idx2 = min(idx1, idx2), max(idx1, idx2)
        counter.add((idx1, idx2))
        data[idx1][idx2] += 1  # type: ignore
    print(f"Heatmap unique combinations: {len(counter)}")

    for i in range(len(ticks)):
        for j in range(i + 1, len(ticks)):
            data[j][i] = None

    result = {
        "tickvals": list(range(len(ticks))),
        "ticktext": ticks,
        "z": data
    }
    return result


def more_eq_than_6_fullreads(item):
    sample, path, data = item
    num_of_full_reads = sum([x["count"] for x in data["nomenclatures"]])
    return num_of_full_reads >= 6


def more_eq_than_10_fullreads(item):
    sample, path, data = item
    num_of_full_reads = sum([x["count"] for x in data["nomenclatures"]])
    return num_of_full_reads >= 10


def more_eq_than_3_allelereads(item):
    sample, path, data = item
    for m in data["modules"]:
        a1_reads = m["allele_1"][4]
        a2_reads = m["allele_2"][4]
        if a1_reads < 3:
            return False
        if a2_reads < 3:
            return False
        a1_pred = m["allele_1"][0]
        a2_pred = m["allele_2"][0]
        # homozygous case
        if a1_pred == a2_pred and a1_reads / 2 < 3:  # covers also a2_reads / 2
            return False

    return True


def more_eq_than_5_allelereads(item):
    sample, path, data = item
    for m in data["modules"]:
        a1_reads = m["allele_1"][4]
        a2_reads = m["allele_2"][4]
        if a1_reads < 5:
            return False
        if a2_reads < 5:
            return False
        a1_pred = m["allele_1"][0]
        a2_pred = m["allele_2"][0]
        # homozygous case
        if a1_pred == a2_pred and a1_reads / 2 < 5:  # covers also a2_reads / 2
            return False

    return True


def generate_df(v):
    data_tmp1 = []
    for sample, path, tmp_data in v:
        data_tmp1.append((
            sample,
            tuple(tmp_data["phased_seqs"]["nomenclature1"]),
            tmp_data["phased_seqs"]["nomenclature1_len"],
            tmp_data["phased_seqs"]["errors1"],
            tuple(tmp_data["phased_seqs"]["nomenclature2"]),
            tmp_data["phased_seqs"]["nomenclature2_len"],
            tmp_data["phased_seqs"]["errors2"]
        ))

    columns = ["sample", "nomenclature1", "nomenclature1_len", "warnings1", "nomenclature2", "nomenclature2_len", "warnings2"]
    df = pd.DataFrame.from_records(data_tmp1, columns=columns)
    df["nomenclature1_join"] = df["nomenclature1"].apply(lambda x: "".join(x))
    df["nomenclature2_join"] = df["nomenclature2"].apply(lambda x: "".join(x))

    # 2025-06-02
    # df["nomenclature1_join"] = df["nomenclature1"].apply(lambda x: x[0])
    # df["nomenclature2_join"] = df["nomenclature2"].apply(lambda x: x[0])
    # df["nomenclature1_len"] = df["nomenclature1_join"].apply(lambda x: len(x))
    # df["nomenclature2_len"] = df["nomenclature2_join"].apply(lambda x: len(x))
    return df


def function2(allele_pairs: list[tuple[int | str, int | str]], max_allele: int) -> dict:
    # z = [
    #     [7   , 0   , 0   , 0   , 0   , 1   , 0   , 0],
    #     [None, 0   , 0   , 0   , 0   , 1   , 0   , 0],
    #     [None, None, 0   , 0   , 0   , 0   , 0   , 0],
    #     [None, None, None, 0   , 35  , 119 , 3   , 3],
    #     [None, None, None, None, 0   , 237 , 3   , 9],
    #     [None, None, None, None, None, 0   , 2   , 0],
    #     [None, None, None, None, None, None, 0   , 1],
    #     [None, None, None, None, None, None, None, 0]
    # ]
    # tickvals = [0, 1, 2, 3, 4, 5, 6, 7]
    # ticktext = ["B", 0, 1, 2, 3, 4, 5, "E"]
    # xlim = 6.5
    z, tickvals, ticktext, xlim = collect_heatmap_data(max_allele, allele_pairs)
    result = {
        "z": z,
        "tickvals": tickvals,
        "ticktext": ticktext,
        "xlim": xlim
    }
    return result


def collect_heatmap_data(max_allele: int, pairs: list[tuple[int | str, int | str]]) -> HeatmapData:
    z: list[list[int]] = [[0] * (max_allele + 1 + 2) for _ in range(max_allele + 1 + 2)]

    for pair in pairs:
        a1 = allele_num(pair[0])
        a2 = allele_num(pair[1])

        if a1 >= 0 and a2 >= 0:
            z[a1 + 1][a2 + 1] += 1
        else:
            if a1 == -1 and a2 == -1:
                z[0][0] += 1
                continue
            if a1 == -1:
                z[0][a2 + 1] += 1
                continue
            if a1 == -3 and a2 == -3:
                z[max_allele + 2][max_allele + 2] += 1
                continue
            if a2 == -3:
                z[a1 + 1][max_allele + 2] += 1
                continue
            raise ValueError(f"{pair} has unexpected value.")

    z_new: list[list[int | None]] = z  # type: ignore
    for i in range(len(z)):
        for j in range(i):
            assert z_new[i][j] == 0, f"z[{i}][{j}] is non-zero"
            z_new[i][j] = None
    tickvals = list(range(max_allele + 3))
    ticktext = ["B"] + list(range(max_allele + 1)) + ["E"]
    xlim = max_allele + 1.5

    heatmap_data = (z_new, tickvals, ticktext, xlim)
    return heatmap_data


def function3(allele_pairs: list[tuple[int | str, int | str]], max_allele: int) -> dict:
    ticktext = ["B"] + list(range(max_allele + 1)) + ["E"]
    tickvals = list(range(len(ticktext)))
    z = [0] * len(ticktext)
    z = collect_alleles(z, allele_pairs, 0)
    z = collect_alleles(z, allele_pairs, 1)

    result = {
        "z": z,
        "tickvals": tickvals,
        "ticktext": ticktext,
    }
    return result


def collect_alleles(hist_data: list[int], allele_pairs: list[tuple[int | str, int | str]], a_idx: int) -> list[int]:
    for pair in allele_pairs:
        allele: int | str = pair[a_idx]
        if allele == 'E':
            hist_data[-1] += 1
            continue
        if allele == 'B':
            hist_data[0] += 1
            continue
        if allele == 'X':
            raise NotImplementedError

        allele = int(allele)
        hist_data[allele + 1] += 1
    return hist_data


def generate_module_data(motif, v, i) -> dict:
    max_allele = MIN_GS
    allele_pairs = []
    rows = []
    for sample, align, data in v:
        a1 = data["modules"][i]["allele_1"][0]
        a2 = data["modules"][i]["allele_2"][0]
        rows.append({
            "sample": sample,
            "a1_pred": a1, "a1_conf": data["modules"][i]["allele_1"][1],
            "a2_pred": a2, "a2_conf": data["modules"][i]["allele_2"][1],
            "conf": data["modules"][i]["stats"][0],
            "spanning_num": data["modules"][i]["reads_spanning"],
            "flanking_num": data["modules"][i]["reads_flanking"],
            "histogram_data": data["modules"][i]["graph_data"][0],
            "alignment_file": align + f"/{motif}.html"
        })
        max_allele = max(max_allele, allele_num(a1))
        max_allele = max(max_allele, allele_num(a2))
        allele_pairs.append((a1, a2))

    heatmap_data = function2(allele_pairs, max_allele)
    histogram_data = function3(allele_pairs, max_allele)

    module = {
        "table": rows,
        "heatmap": heatmap_data,
        "histogram": histogram_data
    }

    return module


def allele_num(x: int | str) -> int:
    if x == "E":
        return -3
    if x == "X":
        return -2
    if x == "B":
        return -1
    return x  # type: ignore


def filter_ticks(ticks: list[str], df: pd.DataFrame, at_least: int) -> list[str]:
    if at_least == 0:
        return ticks
    tmp = create_main_histrogram(df, ticks)
    data = tmp["x"]

    new_ticks = []
    for i, x in enumerate(data):
        if x >= at_least:
            new_ticks.append(ticks[i])

    return new_ticks


def write_structural_nomenclatures(items, motif) -> list[str]:
    result = []
    for item in items:
        sample, path, data = item
        nom1 = " ".join(data["phased_seqs"]["nomenclature1"])
        nom2 = " ".join(data["phased_seqs"]["nomenclature2"])
        line = "\t".join([sample, motif, nom1, nom2])
        result.append(line)
    return result


def main() -> None:
    args = load_arguments()
    motif_dict = collect_jsons(args.inputs, args.output_dir)

    print("Aggregates done. Creating HTMLs.")
    for motif, v in motif_dict.items():
        # 2025-06-04
        # if motif != "DM2" and motif != "HD":
        #     continue

        # 2025-06-02
        # v = list(filter(more_eq_than_6_fullreads, v))
        # v = list(filter(more_eq_than_10_fullreads, v))
        # v = list(filter(more_eq_than_3_allelereads, v))
        # v = list(filter(more_eq_than_5_allelereads, v))

        n_modules = len(v[0][2]["modules"])
        print(f"Creating report of motif {motif} ({n_modules=}). Writting {args.output_dir}/{motif}.html.")

        sequence = convert_modules(v[0][2]["motif_stats"]["modules"])
        df = generate_df(v)
        ticks = create_ticks(df)
        ticks = filter_ticks(ticks, df, args.at_least)

        data: dict[str, Any] = {}
        data["motif_name"] = motif
        data["sequence"] = sequence
        data["main_table"] = json.loads(df.to_json(orient="records"))
        data["main_histogram"] = create_main_histrogram(df, ticks)
        data["main_heatmap"] = create_main_heatmap(df, ticks)

        data["modules"] = [generate_module_data(motif, v, i) for i in range(n_modules)]
        # print(json.dumps(tmp_dict, indent=4))

        template_dir = os.path.dirname(__file__) + "/../templates"
        os.makedirs(f"{args.output_dir}", exist_ok=True)
        env = Environment(loader=FileSystemLoader([template_dir]), trim_blocks=True, lstrip_blocks=True)
        template = env.get_template("population_report.html")
        output = template.render(data=data)
        with open(f"{args.output_dir}/{motif}.html", "w") as f:
            f.write(output)

        output2 = write_structural_nomenclatures(v, motif)
        with open(f"{args.output_dir}/{motif}.structures.tsv", "w") as f:
            f.write("\n".join(output2))

    copy_includes(args.output_dir)


if __name__ == '__main__':
    main()
