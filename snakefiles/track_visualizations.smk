configfile: 'configs/track_visualizations.yaml'

tv_igvjs_bw = config['track_visualizations']['igv_js']['bigwigs']

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
        locus=$(head -n1 {annotation_bed}|cut -f4)
        python scripts/track_list_to_json.py \
        -i {input.track_list} \
        -o {output} \
        -f {input.genome_fa} \
        -g {input.annotation_gtf} \
        -s {annotation_bed} \
        -l $locus
        """
