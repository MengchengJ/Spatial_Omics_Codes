import os
def integrate_report(run_id,report_directory):
    figure_list = [f for f in os.listdir(report_directory) if f.startswith(run_id) and f.endswith('.pdf')]
    with open('seq_report_template.tex',encoding='utf-8') as f:
        template = f.read()
    tex_text = template.replace('RUN_ID',run_id)
    figure_dict = {}
    for figure in figure_list:
        cyc,chn = figure.strip(f'{run_id}_').split('_')[1:3]
        if f'cyc_{cyc}_{chn}' not in figure_dict:
            figure_dict[f'cyc_{cyc}_{chn}'] = \
    f'''
    \\newpage
    \\section{{Cycle {cyc}, {chn.capitalize()}}}
    \\subsection{{亮度分布}}
    \\begin{{figure}}[h!]
        \\begin{{center}}
        \\includegraphics[width=0.65\\textwidth]{{{run_id}_cyc_{cyc}_{chn}_int_dist.pdf}}
        \\end{{center}}
    \\end{{figure}}
    \\subsection{{细胞信号统计}}
    \\begin{{figure}}[h!]
        \\begin{{center}}
        \\includegraphics[width=0.65\\textwidth]{{{run_id}_cyc_{cyc}_{chn}_spot_cell_count.pdf}}
        \\end{{center}}
    \\end{{figure}}
    '''
    figure_text_list = [figure_dict[k] for k in figure_dict]
    tex_text = tex_text.replace('CONTEXT','\n'.join(figure_text_list))
    file_name = f'{run_id}_report.tex'
    with open(os.path.join(report_directory,file_name),'w',encoding='utf-8') as f:
        f.write(tex_text)
    os.chdir(report_directory)
    for _ in range(2):
        os.system(f'xelatex -interaction=batchmode -halt-on-error {file_name}')

if __name__ == "__main__":
    pass