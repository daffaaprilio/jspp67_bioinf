# workflows/kegg_scraping/kegg_ftp.smk

configfile: "./config.yaml"

KEGG_USER   = config['user']
KEGG_PASS   = config['pass']
# Define URLs for KEGG FTP downloads
KO_URL          = config['url_ko']
SBI_URL         = config['url_sbi']
MAP_URL         = config['url_pathway']

"""
Download & prepare the following from KEGG FTP:
    ko.tar.gz (EC to KO mapping)
    sbi (KO to organism mapping)
    map (pathway to KO list)
"""

rule all:
    input:
        '../../data/reference/KEGG/ko.tar.gz',
        '../../data/reference/KEGG/sbi',
        '../../data/reference/KEGG/map.tar.gz'

rule obtain_ko:
    params:
        user        = KEGG_USER,
        password    = KEGG_PASS,
        url         = KO_URL
    output:
        archive     = '../../data/reference/KEGG/ko.tar.gz',
        outdir      = '../../data/reference/KEGG/ko'
    shell:
        """
        wget --user {params.user} --password {params.password} -O {output.archive} {params.url}
        mkdir -p {output.outdir}
        """

rule obtain_sbi:
    params:
        user        = KEGG_USER,
        password    = KEGG_PASS,
        url         = SBI_URL
    output:
        outdir      = directory('../../data/reference/KEGG/sbi'),
    shell:
        """
        mkdir -p {output.outdir}
        wget -m -np -e robots=off \
            --user {params.user} --password {params.password} \
            --no-host-directories \
            --cut-dirs=4 \
            -P {output.outdir} \
            {params.url}
        """

rule obtain_map:
    params:
        user        = KEGG_USER,
        password    = KEGG_PASS,
        url         = MAP_URL
    output:
        archive     = '../../data/reference/KEGG/map.tar.gz',
        outdir      = '../../data/reference/KEGG/map'
    shell:
        """
        wget --user {params.user} --password {params.password} -O {output.archive} {params.url}
        mkdir -p {output.outdir}
        """
