configfile: 'configs/track_visualizations.yaml'

tv_igvjs_bw = config['track_visualizations']['igv_js']['bigwigs']
tv_igvjs_redirect_template = config['track_visualizations']['igv_js']['session_redirect_template']
tv_igvjs_redirect_string = config['track_visualizations']['igv_js']['session_redirect_string']

rule track_visualizations_igv_js_bigwigs:
    input:
        track_list = tv_igvjs_bw['input'],
        genome_fa = prefix+'genome.fa',
        genome_index = prefix+'genome.fa.fai',
        annotation_gtf = prefix+'annotation.gtf.gz',
        annotation_index = prefix+'annotation.gtf.gz.tbi',
        annotation_bed = prefix+'annotation.simplified.bed'
    output:
        tv_igvjs_bw['output']
    conda:
        '../envs/seaborn.yaml'
    shell:
        """
        
        mkdir -p $(dirname {output})
        locus=$(head -n1 {input.annotation_bed}|cut -f4)
        python scripts/track_list_to_json.py \
        -i {input.track_list} \
        -o {output} \
        -f {input.genome_fa} \
        -g {input.annotation_gtf} \
        -s {input.annotation_bed} \
        -l $locus
        """

rule track_visualizations_igv_js_to_html:
    input:
        json = '{filepath}.igv.json',
        html = tv_igvjs_redirect_template
    output:
        '{filepath}.igv.html'
    params:
        session_redirect_string = tv_igvjs_redirect_string
    run:
        with open(input.html) as template:
            output_html = template.read().replace(params.session_redirect_string, input.json)
            with open(output[0],'w') as outfile:
                outfile.write(output_html)
