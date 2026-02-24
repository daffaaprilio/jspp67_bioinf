# workflows/kegg_scraping/kegg_ftp.smk

## next update: add unzipping (tar gz) to the downloaded files

"""
Download & prepare the following from KEGG FTP:
    ko.tar.gz (EC to KO mapping)
    KEGG_SP (KO to organism mapping)
    map (pathway to KO list)
"""

configfile: "./config-kegg_ftp.yaml"

KEGG_USER           = config['user']
KEGG_PASS           = config['pass']
# Define URLs for KEGG FTP downloads
KO_URL              = config['url_ko']
KEGG_SP_BASE_URL    = config['url_kegg_sp_base']
MAP_URL             = config['url_pathway']

ORGANISMS = config['organisms']

rule all:
    input:
        '../../data/reference/KEGG/ko.tar.gz',
        expand('../../data/reference/KEGG/{organism}', organism=ORGANISMS),
        '../../data/reference/KEGG/map.tar.gz'

rule obtain_ko:
    params:
        user        = KEGG_USER,
        password    = KEGG_PASS,
        url         = KO_URL
    output:
        archive     = '../../data/reference/KEGG/ko.tar.gz'
    shell:
        """
        wget --user {params.user} --password {params.password} -O {output.archive} {params.url}
        """

rule obtain_kegg_sp:
    params:
        user        = KEGG_USER,
        password    = KEGG_PASS,
        url         = lambda wc: f"{KEGG_SP_BASE_URL}{wc.organism}/"
    output:
        outdir      = directory('../../data/reference/KEGG/{organism}'),
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
        archive     = '../../data/reference/KEGG/map.tar.gz'
    shell:
        """
        wget --user {params.user} --password {params.password} -O {output.archive} {params.url}
        """
