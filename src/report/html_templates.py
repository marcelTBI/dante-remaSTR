import gzip
import numpy as np
import re

import pandas as pd

from src.postfilter import PostFilter

contents = """
<table class="tg" id="content-tg">
    <thead>
        <tr>
            <th class="tg-s6z2">Motif</th>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>
"""

make_datatable_string = """
<script>
    $(document).ready( function () {{
    $('#content-tg').DataTable();
}} );
</script>
"""

content_string = """ <tr>
    <td class="tg-s6z2">
        <a href="#{motif_name}" class="tablinks" onclick="openTab(event, '{motif_name}'); $('.{motif_name}').trigger('content-change');">{motif}</a>
    </td>
</tr>"""

content_string_empty = ""

row_string = """  <tr>
    <td class="tg-s6z2">{motif_name}</td>
    <td class="tg-s6z2">{motif_nomenclature}</td>
    <td class="tg-s6z2">{allele1}</td>
    <td class="tg-s6z2">{conf_allele1}</td>
    <td class="tg-s6z2">{reads_a1}</td>
    <td class="tg-s6z2">{indels_a1}</td>
    <td class="tg-s6z2">{mismatches_a1}</td>
    <td class="tg-s6z2">{allele2}</td>
    <td class="tg-s6z2">{conf_allele2}</td>
    <td class="tg-s6z2">{reads_a2}</td>
    <td class="tg-s6z2">{indels_a2}</td>
    <td class="tg-s6z2">{mismatches_a2}</td>
    <td class="tg-s6z2">{confidence}</td>
    <td class="tg-s6z2">{indels}</td>
    <td class="tg-s6z2">{mismatches}</td>
    <td class="tg-s6z2">{quality_reads}</td>
    <td class="tg-s6z2">{one_primer_reads}</td>
  </tr>"""

nomenclature_string = """
<tr>
    <td>{count}</td>
    <td>{ref}</td>
    {parts}
</tr>
"""

motif_summary = """
<div class="tabcontent" id="{motif_id}" style="display: none">
<h2 class="summary_nomenclatures">Nomenclatures</h2>
<table class="nomtg">
    <tbody>
        {nomenclatures}
    </tbody>
</table>
<h2 class="summary">Summary table</h2>
<table class="tg" id="tg-{motif_id}">
    <thead>
        <tr>
            <th class="tg-s6z2" rowspan="3">Motif</th>
            <th class="tg-s6z2" rowspan="3">Sequence</th>
            <th class="tg-s6z2" colspan="5">Allele 1</th>
            <th class="tg-s6z2" colspan="5">Allele 2</th>
            <th class="tg-s6z2" rowspan="3">Overall<br>confidence</th>
            <th class="tg-s6z2" colspan="2">Errors per read</th>
            <th class="tg-s6z2" colspan="2">Reads</th>
        </tr>
        <tr>
            <td class="tg-smaller" rowspan="2">prediction</td>
            <td class="tg-smaller" rowspan="2">confidence</td>
            <td class="tg-smaller" rowspan="2">reads</td>
            <td class="tg-smaller" colspan="2">Errors per read</td>
            <td class="tg-smaller" rowspan="2">prediction</td>
            <td class="tg-smaller" rowspan="2">confidence</td>
            <td class="tg-smaller" rowspan="2">reads</td>
            <td class="tg-smaller" colspan="2">Errors per read</td>
            <td class="tg-smaller" rowspan="2">indel</td>
            <td class="tg-smaller" rowspan="2">mismatch</td>
            <td class="tg-smaller" rowspan="2">full</td>
            <td class="tg-smaller" rowspan="2">partial</td>
        </tr>
        <tr>
            <td class="tg-smaller">indel</td>
            <td class="tg-smaller">mismatch</td>
            <td class="tg-smaller">indel</td>
            <td class="tg-smaller">mismatch</td>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>

<script>
    $(document).ready( function () {{
    $('#tg-{motif_id}').DataTable({{scrollX: true}});
}} );
</script>
<p><a href="#content">Back to content</a></p>
{motifs}
</div>
"""

motif_summary_static = """
<div style="margin: 25px 0 25px 0">
<h2 class="summary">Summary table</h2>
<table class="tg">
    <thead>
        <tr>
            <th class="tg-s6z2" rowspan="3">Motif</th>
            <th class="tg-s6z2" rowspan="3">Sequence</th>
            <th class="tg-s6z2" colspan="5">Allele 1</th>
            <th class="tg-s6z2" colspan="5">Allele 2</th>
            <th class="tg-s6z2" rowspan="3">Overall<br>confidence</th>
            <th class="tg-s6z2" colspan="2">Errors per read</th>
            <th class="tg-s6z2" colspan="2">Reads</th>

        </tr>
        <tr>
            <td class="tg-s6z2" rowspan="2">prediction</td>
            <td class="tg-s6z2" rowspan="2">confidence</td>
            <td class="tg-s6z2" rowspan="2">reads</td>
            <td class="tg-s6z2" colspan="2">Errors per read</td>
            <td class="tg-s6z2" rowspan="2">prediction</td>
            <td class="tg-s6z2" rowspan="2">confidence</td>
            <td class="tg-s6z2" rowspan="2">reads</td>
            <td class="tg-s6z2" colspan="2">Errors per read</td>
            <td class="tg-s6z2" rowspan="2">indel</td>
            <td class="tg-s6z2" rowspan="2">mismatch</td>
            <td class="tg-s6z2" rowspan="2">full</td>
            <td class="tg-s6z2" rowspan="2">partial</td>
        </tr>
        <tr>
            <td class="tg-s6z2">indel</td>
            <td class="tg-s6z2">mismatch</td>
            <td class="tg-s6z2">indel</td>
            <td class="tg-s6z2">mismatch</td>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>
</div>
"""

motif_stringb64 = """
<h2 id="data-{motif_name}">{motif}</h2>
<p>{sequence}</p>
Nomenclatures:<br>
<table class="nomtg">
    <tbody>
        {nomenclatures}
    </tbody>
</table>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist pic100 {motif_id}" id="hist-{motif_name}"></div>
        </td>
        <td colspan="1">
            <div class="pcol pic100 {motif_id}" id="pcol-{motif_name}"></div>
        </td>
    </tr>
</table>

<script>
    {{
        let hist_data = {motif_reps};
        let pcol_data = {motif_pcolor};
        
        let updateGraph = () => {{
            if (document.getElementById('{motif_id}').style.display === 'block') {{
                hist_data['layout'] = {{...hist_data['layout'], width: (window.innerWidth-50) * 0.6, height: (window.innerWidth-50) * 0.35,
                    legend: {{...hist_data['layout']['legend'], x: 0.85, y: 0.85 }}}};
                pcol_data['layout'] = {{...pcol_data['layout'], width: (window.innerWidth-50) * 0.4, height: (window.innerWidth-50) * 0.35}};
                Plotly.react('hist-{motif_name}', hist_data);
                Plotly.react('pcol-{motif_name}', pcol_data);
            }}
        }};
        
        $(document).ready(function() {{
            $('.{motif_id}').bind("content-change", updateGraph);
            window.addEventListener('resize', updateGraph, true);
        }});
    }}
</script>

<p><a href="{alignment}">Link to alignments</a></p>
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_reponly = """
<h2 id="data-{motif_name}">{motif}</h2>
<p>{sequence}</p>
Nomenclatures:<br>
<table class="nomtg">
    <tbody>
        {nomenclatures}
    </tbody>
</table>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist {motif_id}" id="hist2d-{motif_name}"></div>
        </td>
    </tr>
</table>

<script>
    {{
        let hist_data = {motif_reps};
        
        let updateGraph = () => {{
            if (document.getElementById('{motif_id}').style.display === 'block') {{
                hist_data['layout'] = {{...hist_data['layout'], width: (window.innerWidth-50) * 0.5, height: (window.innerWidth-50) * 0.45}};
                Plotly.react('hist2d-{motif_name}', hist_data);
            }}
        }};
        
        $(document).ready(function() {{
            $('.{motif_id}').bind("content-change", updateGraph);
            window.addEventListener('resize', updateGraph, true);
        }});
    }}
</script>

<p><a href="{alignment}">Link to alignments</a></p>
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_static = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
Nomenclatures:<br>
<table class="nomtg">
    <tbody>
        {nomenclatures}
    </tbody>
</table>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist pic100 {motif_id}" id="hist-{motif_name}"></div>
        </td>
        <td colspan="1">
            <div class="pcol pic100 {motif_id}" id="pcol-{motif_name}"></div>
        </td>
    </tr>
</table>

<script>
    {{
        let hist_data = {motif_reps};
        let pcol_data = {motif_pcolor};
        
        let updateGraph = () => {{
            hist_data['layout'] = {{...hist_data['layout'], width: (window.innerWidth-50) * 0.6, height: (window.innerWidth-50) * 0.35,
                legend: {{...hist_data['layout']['legend'], x: 0.85, y: 0.85 }}}};
            pcol_data['layout'] = {{...pcol_data['layout'], width: (window.innerWidth-50) * 0.4, height: (window.innerWidth-50) * 0.35}};
            Plotly.react('hist-{motif_name}', hist_data);
            Plotly.react('pcol-{motif_name}', pcol_data);
        }};
        updateGraph();
        $(document).ready(function() {{
            window.addEventListener('resize', updateGraph, true);
        }});
    }}
</script>

<p><a href="{alignment}">Link to alignments</a></p>
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_reponly_static = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
Nomenclatures:<br>
<table class="nomtg">
    <tbody>
        {nomenclatures}
    </tbody>
</table>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist pic50 {motif_id}" id="hist2d-{motif_name}"></div>
        </td>
    </tr>
</table>

<script>
    {{
        let hist_data = {motif_reps};
        
        let updateGraph = () => {{
            hist_data['layout'] = {{...hist_data['layout'], width: (window.innerWidth-50) * 0.5, height: (window.innerWidth-50) * 0.45}};
            Plotly.react('hist2d-{motif_name}', hist_data);
        }};
        updateGraph();
        $(document).ready(function() {{
            window.addEventListener('resize', updateGraph, true);
        }});
    }}
</script>

<p><a href="#content">Back to content</a></p>
"""

motif_string_empty_old = """<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
NOT AVAILABLE
<p><a href="#content">Back to content</a></p>
"""

motif_string_empty = ""

alignment_string = """
  <p>{sequence}</p>
  {alignment}
  <hr>
"""

align_vis = """
  <details>
    <summary>{display_text}</summary>
    <div id="A{name}" class="align">press "Run with JS"</div>
    <script>
        var fasta = `{fasta}`;
        var seqs = msa.io.fasta.parse(fasta);
        var opts = {{
            el: document.getElementById("A{name}"),
            vis: {{
                conserv: false,
                metaIdentity: true,
                overviewbox: true,
                seqlogo: {seq_logo}
            }},
            seqs: seqs,
            colorscheme: {{"scheme": "nucleotide"}},
            // smaller menu for JSBin
            menu: "small",
            bootstrapMenu: true
        }};
        var m = new msa.msa(opts);
        m.render()
    </script>
  </details>
"""


def highlight_subpart(seq: str, highlight: int | list[int]) -> tuple[str, str]:
    """
    Highlights subpart of a motif sequence
    :param seq: str - motif sequence
    :param highlight: int/list(int) - part ot highlight
    :return: str, str - motif sequence with highlighted subpart, highlighted subpart
    """
    if highlight is None:
        return seq, ''

    str_part = []
    highlight = np.array(highlight)
    split = [f'{s}]' for s in seq.split(']') if s != '']
    for h in highlight:
        str_part.append(split[h])
        split[h] = f'<b><u>{split[h]}</u></b>'
    return ''.join(split), ''.join(str_part)


def float_to_str(c: float | str, percents: bool = False, decimals: int = 1) -> str:
    """
    Convert float confidence to string.
    :param c: float/str - confidence
    :param percents: bool - whether to output as a percents or not
    :param decimals: int - how many decimals to round to
    :return: str - converted to string
    """
    if isinstance(c, float):
        return f'{c * 100: .{decimals}f}%' if percents else f'{c: .{decimals}f}'
    return c


def generate_row(sequence: str, result: pd.Series, postfilter: PostFilter) -> str:
    """
    Generate rows of a summary table in html report.
    :param sequence: str - motif sequence
    :param result: pd.Series - result row to convert to table
    :param postfilter: PostFilter - postfilter dict from config
    :return: str - html string with rows of the summary table
    """
    highlight = list(map(int, str(result['repetition_index']).split('_')))
    sequence, subpart = highlight_subpart(sequence, highlight)

    # shorten sequence:
    keep = 10
    first = sequence.find(',')
    last = sequence.rfind(',')
    smaller_seq = sequence if first == -1 else '...' + sequence[first - keep:last + keep + 1] + '...'

    # errors:
    errors = f'{postfilter.max_rel_error * 100:.0f}%'
    if postfilter.max_abs_error is not None:
        errors += f' (abs={postfilter.max_abs_error})'

    # fill templates:
    updated_result = {'conf_allele1': float_to_str(result['conf_allele1'], percents=True),
                      'conf_allele2': float_to_str(result['conf_allele2'], percents=True),
                      'confidence': float_to_str(result['confidence'], percents=True), 'motif_nomenclature': smaller_seq,
                      'indels': float_to_str(result['indels'], decimals=2), 'mismatches': float_to_str(result['mismatches'], decimals=2),
                      'indels_a1': float_to_str(result['indels_a1'], decimals=2), 'mismatches_a1': float_to_str(result['mismatches_a1'], decimals=2),
                      'indels_a2': float_to_str(result['indels_a2'], decimals=2), 'mismatches_a2': float_to_str(result['mismatches_a2'], decimals=2)}
    return row_string.format(**{**result, **updated_result})


def get_alignment_name(alignment_file: str, allele: int) -> str:
    """
    Get alignment file of subpart of the alignment with allele count specified.
    :param alignment_file: str - alignment file name
    :param allele: int - allele repetition count
    :return: str - alignment file name for the subpart of alignment
    """
    # find where is .fasta
    fasta_index = alignment_file.rfind('.fasta')
    # insert '_aX' before .fasta
    return alignment_file[:fasta_index] + '_a' + str(allele) + alignment_file[fasta_index:]


def generate_motifb64(sequence: str, result: pd.Series, repetition: str, pcolor: str | None, alignment: str | None, filtered_alignment: str | None,
                      filtered_left_alignment: str | None, filtered_right_alignment: str | None, nomenclature_lines: list[str],
                      postfilter: PostFilter, static: bool = False) -> tuple[str, str, tuple[str, str]]:
    """
    Generate part of a html report for each motif.
    :param sequence: str - motif sequence
    :param repetition: str - filename of repetitions figures
    :param pcolor: str - filename of pcolor figures
    :param alignment: str/None - filename of alignment file
    :param filtered_alignment: str/None - filename of filtered alignment file
    :param filtered_left_alignment: str/None - filename of filtered left alignment file
    :param filtered_right_alignment: str/None - filename of filtered right alignment file
    :param nomenclature_lines: list[str] - list of nomenclature strings
    :param postfilter: PostFilter - postfilter arguments
    :param static: bool - generate static code?
    :return: (str, str) - content and main part of the html report for motifs
    """
    # prepare and generate alignments
    highlight = list(map(int, str(result['repetition_index']).split('_')))
    sequence, subpart = highlight_subpart(sequence, highlight)
    motif = result['motif_name']
    motif_name = f'{result["motif_name"].replace("/", "_")}_{",".join(map(str, highlight)) if highlight is not None else "mot"}'
    motif_clean = re.sub(r'[^\w_]', '', motif_name)
    motif_clean_id = motif_clean.rsplit('_', 1)[0] if highlight == [1] else motif_clean  # trick to solve static html
    align_html_a1 = ''
    align_html_a2 = ''

    a1 = result['allele1']
    a2 = result['allele2']
    c = result['confidence']
    c1 = result['conf_allele1']
    c2 = result['conf_allele2']
    if (a1 == 'B' and a2 == 'B') or (a1 == 0 and a2 == 0):
        result = f'BG {float_to_str(c, percents=True)}'
    else:
        result = f'{str(a1):2s} ({float_to_str(c1, percents=True)}) {str(a2):2s} ({float_to_str(c2, percents=True)}) total {float_to_str(c, percents=True)}'
        if alignment is not None:
            align_html_a1 = generate_alignment(f'{motif_clean}_{str(a1)}', get_alignment_name(alignment, a1), motif_clean.split('_')[0],
                                               f'Allele 1 ({str(a1):2s}) alignment visualization')
            if a1 != a2:
                align_html_a2 = generate_alignment(f'{motif_clean}_{str(a2)}', get_alignment_name(alignment, a2), motif_clean.split('_')[0],
                                                   f'Allele 2 ({str(a2):2s}) alignment visualization')

    # errors:
    errors = f'{postfilter.max_rel_error * 100:.0f}%'
    if postfilter.max_abs_error is not None:
        errors += f' (abs={postfilter.max_abs_error})'

    # return content and picture parts:
    motif_templates = {'static': {'pcol': motif_stringb64_static, 'no-pcol': motif_stringb64_reponly_static},
                       'dynamic': {'pcol': motif_stringb64, 'no-pcol': motif_stringb64_reponly}}

    if repetition is not None:
        reps = open(repetition, 'r').read()
        align_html = generate_alignment(motif_clean, alignment, motif_clean.split('_')[0], 'Spanning reads alignment visualization')
        filt_align_html = generate_alignment(motif_clean + '_filtered', filtered_alignment, motif_clean.split('_')[0],
                                             'Partial reads alignment visualization', seq_logo=False)
        left_align_html = generate_alignment(motif_clean + '_filtered_left', filtered_left_alignment, motif_clean.split('_')[0],
                                             'Left flank reads alignment visualization')
        right_align_html = generate_alignment(motif_clean + '_filtered_right', filtered_right_alignment, motif_clean.split('_')[0],
                                              'Right flank reads alignment visualization')

        # select template
        motif_template = motif_templates['static' if static else 'dynamic']['no-pcol' if pcolor is None else 'pcol']

        # read pcolor if available
        pcol = '' if pcolor is None else open(pcolor, 'r').read()

        # return filled valid template
        return (content_string.format(motif_name=motif_clean.rsplit('_', 1)[0], motif=motif),
                motif_template.format(post_bases=postfilter.min_rep_len, post_reps=postfilter.min_rep_cnt, motif_name=motif_clean_id,
                                      motif_id=motif_clean.rsplit('_', 1)[0], motif=motif, motif_reps=reps, result=result, motif_pcolor=pcol,
                                      alignment=f'{motif.replace(" ", "%20")}/alignments.html', sequence=sequence, errors=errors,
                                      nomenclatures='\n'.join(nomenclature_lines)),
                (motif, alignment_string.format(sequence=sequence, alignment=align_html + align_html_a1 + align_html_a2 +
                                                                             filt_align_html + left_align_html + right_align_html)))
    else:
        return (content_string_empty.format(motif_name=motif_clean.rsplit('_', 1)[0], motif=motif),
                motif_string_empty.format(post_bases=postfilter.min_rep_len, post_reps=postfilter.min_rep_cnt,
                                          motif_name=motif_clean, motif=motif, sequence=sequence, errors=errors),
                (motif, ''))


def generate_alignment(motif: str, alignment_file: str, motif_id: str, display_text: str = 'Click to toggle alignment visualization',
                       seq_logo: bool = True) -> str:
    """
    Generate HTML code for the fancy alignment.
    :param motif: str - name of the motif
    :param alignment_file: str - filename of the alignment file
    :param motif_id: str - motif identification
    :param display_text: str - string to display when the alignment is hidden
    :param seq_logo: bool - display sequence logo?
    :return: str - code of the fancy alignment
    """
    if alignment_file is None:
        return ''

    try:
        with gzip.open(alignment_file, 'rt') if alignment_file.endswith('.gz') else open(alignment_file) as f:
            string = f.read()
        string = string[:string.find('#')]

        return align_vis.format(fasta=string, name=motif, motif_id=motif_id, display_text=display_text, seq_logo='true' if seq_logo else 'false')
    except (IOError, TypeError, AttributeError):
        return ''
