#!/usr/bin/env python3
"""
Methylation Landscape EDA — standalone script
==============================================
Headless equivalent of methylation_landscape.ipynb.
Produces seven PDF figures in --out-dir.

Usage (from project root, no args required):
    python workflows/methylation_landscape/run_eda.py

Usage with explicit paths:
    python workflows/methylation_landscape/run_eda.py \
        --wdir /home/daffa/Work/2026/02-JSPP67 \
        --out-dir results/methylation_landscape

Background:
    nohup python workflows/methylation_landscape/run_eda.py \
        > results/methylation_landscape/run_eda.log 2>&1 &
"""

import argparse
import logging
import subprocess
import sys
import tempfile
from functools import reduce
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use('Agg')                      # must be before pyplot import
import matplotlib.pyplot as plt            # noqa: E402
import matplotlib.ticker as mticker        # noqa: E402
import numpy as np                         # noqa: E402
import pandas as pd                        # noqa: E402
import pyfaidx                             # noqa: E402
import seaborn as sns                      # noqa: E402
from sklearn.decomposition import PCA      # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
log = logging.getLogger(__name__)

# ── Column names (modkit pileup bedMethyl, 18 columns) ───────────────────────
BEDMETHYL_COLS = [
    'chrom', 'start', 'end', 'mod_code', 'score', 'strand',
    'start2', 'end2', 'color',
    'N_valid_cov', 'pct_mod', 'N_mod',
    'N_canonical', 'N_other_mod', 'N_del', 'N_fail', 'N_diff', 'N_nocall',
]

SAMPLE_META = {
    'SBC4':  {'label': 'SBC4\n(TAA med)',   'color': '#4CAF50'},
    'SBC10': {'label': 'SBC10\n(TAA high)', 'color': '#2196F3'},
    'SBC11': {'label': 'SBC11\n(TAA low)',  'color': '#F44336'},
    'SBC23': {'label': 'SBC23\n(TAA med)',  'color': '#8BC34A'},
}


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_bedmethyl(path, usecols=None, extra_dtypes=None):
    dtypes = {'chrom': str, 'start': int, 'end': int, 'strand': str,
              'N_valid_cov': int, 'pct_mod': float, 'N_mod': int}
    if extra_dtypes:
        dtypes.update(extra_dtypes)
    return pd.read_csv(
        path, sep='\t', header=None, names=BEDMETHYL_COLS,
        usecols=usecols, dtype=dtypes,
    )


def discover_available(results_dir, samples, contexts):
    """Return dict {context: [samples whose bedMethyl exists]}."""
    avail = {}
    for ctx in contexts:
        avail[ctx] = [
            s for s in samples
            if (results_dir / ctx / s / f'{s}.bedMethyl').exists()
        ]
    return avail


def classify_contexts(df, fasta):
    """
    Add 'context' (CG / CHG / CHH) to a bedMethyl DataFrame.
    df must have a default RangeIndex and columns: chrom, start, strand.
    """
    ctx_arr   = np.full(len(df), 'CHH', dtype='U3')
    avail_chr = set(fasta.keys())

    for chrom, grp in df.groupby('chrom', sort=False):
        if chrom not in avail_chr:
            continue
        seq  = np.frombuffer(
            fasta[chrom][:].seq.upper().encode('ascii'), dtype='S1')
        clen = len(seq)

        plus_idx = grp.index[grp['strand'] == '+'].values
        if len(plus_idx):
            p  = df.loc[plus_idx, 'start'].values
            n1 = seq[np.clip(p + 1, 0, clen - 1)]
            n2 = seq[np.clip(p + 2, 0, clen - 1)]
            ctx = np.where(n1 == b'G', 'CG', np.where(n2 == b'G', 'CHG', 'CHH'))
            ctx[p + 2 >= clen] = 'CHH'
            ctx_arr[plus_idx]  = ctx

        minus_idx = grp.index[grp['strand'] == '-'].values
        if len(minus_idx):
            p  = df.loc[minus_idx, 'start'].values
            m1 = seq[np.clip(p - 1, 0, clen - 1)]
            m2 = seq[np.clip(p - 2, 0, clen - 1)]
            ctx = np.where(m1 == b'C', 'CG', np.where(m2 == b'C', 'CHG', 'CHH'))
            ctx[p < 2] = 'CHH'
            ctx_arr[minus_idx] = ctx

    return df.assign(context=ctx_arr)


# ─────────────────────────────────────────────────────────────────────────────
# Section 1 — Average methylation per context (CG / CHG / CHH)
# ─────────────────────────────────────────────────────────────────────────────

def section1(results_dir, available, ref_fasta, out_dir, sample_order):
    log.info('=== Section 1: Average methylation by context ===')
    samples = available.get('allC', [])
    if not samples:
        log.warning('No allC bedMethyl files available — skipping Section 1.')
        return {}

    # Load allC data
    allC_data = {}
    for sample in samples:
        fpath = results_dir / 'allC' / sample / f'{sample}.bedMethyl'
        log.info(f'  Loading {sample} ...')
        df = load_bedmethyl(fpath)
        allC_data[sample] = df[['chrom', 'start', 'end', 'strand',
                                  'N_valid_cov', 'pct_mod', 'N_mod']].copy()
        log.info(f'    {len(df):,} rows | {df["chrom"].nunique()} chromosomes')

    # Context classification
    log.info(f'  Opening reference FASTA: {ref_fasta}')
    genome = pyfaidx.Fasta(str(ref_fasta), build_index=False)
    allC_ctx = {}
    for sample, df in allC_data.items():
        log.info(f'  Classifying contexts for {sample} ({len(df):,} positions) ...')
        df_ri = df.reset_index(drop=True)
        allC_ctx[sample] = classify_contexts(df_ri, genome)
        vc = allC_ctx[sample]['context'].value_counts()
        log.info(f'    CG: {vc.get("CG",0):,}  CHG: {vc.get("CHG",0):,}  CHH: {vc.get("CHH",0):,}')

    # Weighted mean methylation per context per sample
    records = []
    for sample, df in allC_ctx.items():
        for ctx, grp in df.groupby('context'):
            n_cov = grp['N_valid_cov'].sum()
            n_mod = grp['N_mod'].sum()
            wmean = (n_mod / n_cov * 100) if n_cov > 0 else 0.0
            records.append({'sample': sample, 'context': ctx,
                             'mean_meth_pct': wmean, 'n_sites': len(grp)})
    ctx_summary = pd.DataFrame(records)

    # Plot
    CTX_ORDER  = ['CG', 'CHG', 'CHH']
    CTX_COLORS = {'CG': '#2196F3', 'CHG': '#FF9800', 'CHH': '#4CAF50'}
    sorder = [s for s in sample_order if s in allC_ctx]

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(sorder))
    w = 0.25
    for i, ctx in enumerate(CTX_ORDER):
        vals = []
        for s in sorder:
            row = ctx_summary[
                (ctx_summary['sample'] == s) & (ctx_summary['context'] == ctx)]
            vals.append(row['mean_meth_pct'].values[0] if len(row) else 0.0)
        ax.bar(x + (i - 1) * w, vals, w, label=ctx, color=CTX_COLORS[ctx],
               edgecolor='white', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([SAMPLE_META.get(s, {}).get('label', s) for s in sorder])
    ax.set_ylabel('Weighted mean methylation (%)')
    ax.set_title('Average methylation by sequence context (allC bedMethyl)')
    ax.legend(title='Context', frameon=False)
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())
    sns.despine(ax=ax)
    plt.tight_layout()
    out = out_dir / 'fig1_avg_methylation_by_context.pdf'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log.info(f'  Saved: {out}')

    return allC_data, allC_ctx


# ─────────────────────────────────────────────────────────────────────────────
# Section 2 — Methylation at genomic features
# ─────────────────────────────────────────────────────────────────────────────

def section2(results_dir, available, ref_gff, out_dir, tmpdir, sample_order,
             promoter_bp=2000):
    log.info('=== Section 2: Methylation at genomic features ===')
    samples = available.get('CpG', [])
    if not samples:
        log.warning('No CpG bedMethyl files available — skipping Section 2.')
        return {}

    # Parse GFF3
    GFF_COLS = ['seqid', 'source', 'type', 'start', 'end',
                'score', 'strand', 'phase', 'attributes']
    gff   = pd.read_csv(ref_gff, sep='\t', comment='#', header=None,
                        names=GFF_COLS, dtype={'start': int, 'end': int})
    genes = gff[gff['type'] == 'gene'].copy()
    log.info(f'  {len(genes):,} gene features loaded from GFF3')

    gene_bed = pd.DataFrame({
        'chrom':  genes['seqid'].values,
        'start':  genes['start'].values - 1,
        'end':    genes['end'].values,
        'name':   'gene_body', 'score': 0,
        'strand': genes['strand'].values,
    })

    plus_mask   = genes['strand'].values == '+'
    tss         = np.where(plus_mask, genes['start'].values - 1, genes['end'].values)
    promo_start = np.where(plus_mask, np.maximum(0, tss - promoter_bp), tss)
    promo_end   = np.where(plus_mask, tss, tss + promoter_bp)
    promoter_bed = pd.DataFrame({
        'chrom': genes['seqid'].values, 'start': promo_start, 'end': promo_end,
        'name': 'promoter', 'score': 0, 'strand': genes['strand'].values,
    })
    promoter_bed = promoter_bed[promoter_bed['start'] < promoter_bed['end']].copy()

    # Chromosome name check
    _peek = pd.read_csv(
        results_dir / 'CpG' / samples[0] / f'{samples[0]}.bedMethyl',
        sep='\t', header=None, usecols=[0], nrows=2000, names=['chrom'],
    )
    shared = set(_peek['chrom'].unique()) & set(gene_bed['chrom'].unique())
    if not shared:
        log.warning('No shared chromosomes between bedMethyl and GFF3! '
                    'Section 2 will produce empty feature intersections.')
    else:
        log.info(f'  Chromosome names compatible ({len(shared)} shared)')

    gene_bed_path     = tmpdir / 'genes.bed'
    promoter_bed_path = tmpdir / 'promoters.bed'
    gene_bed.sort_values(['chrom', 'start']).to_csv(
        gene_bed_path, sep='\t', header=False, index=False)
    promoter_bed.sort_values(['chrom', 'start']).to_csv(
        promoter_bed_path, sep='\t', header=False, index=False)

    def intersect_row_indices(a_bed, b_bed):
        res = subprocess.run(
            ['bedtools', 'intersect', '-a', str(a_bed), '-b', str(b_bed), '-u'],
            capture_output=True, text=True, check=True,
        )
        if not res.stdout.strip():
            return set()
        out = pd.read_csv(StringIO(res.stdout), sep='\t', header=None)
        return set(out[3].astype(int))

    feat_data = {}
    for sample in samples:
        log.info(f'  Assigning feature types for {sample} ...')
        fpath = results_dir / 'CpG' / sample / f'{sample}.bedMethyl'
        bm = load_bedmethyl(fpath, usecols=[0, 1, 2, 9, 10])
        bm = bm[['chrom', 'start', 'end', 'N_valid_cov', 'pct_mod']]
        bm = bm.sort_values(['chrom', 'start']).reset_index(drop=True)
        bm['row_idx'] = bm.index
        tmp_bed = tmpdir / f'{sample}_cpg.tmp.bed'
        bm[['chrom', 'start', 'end', 'row_idx']].to_csv(
            tmp_bed, sep='\t', header=False, index=False)

        gene_idx  = intersect_row_indices(tmp_bed, gene_bed_path)
        promo_idx = intersect_row_indices(tmp_bed, promoter_bed_path)

        bm['feature_type'] = 'intergenic'
        bm.loc[bm['row_idx'].isin(gene_idx),  'feature_type'] = 'gene_body'
        bm.loc[bm['row_idx'].isin(promo_idx), 'feature_type'] = 'promoter'
        feat_data[sample] = bm.drop(columns='row_idx')

        vc = feat_data[sample]['feature_type'].value_counts()
        total = vc.sum()
        log.info('    ' + '  '.join(
            f'{k}: {v:,} ({v/total*100:.1f}%)' for k, v in vc.items()))

    # Violin plot
    feat_plot_df = pd.concat(
        [df.assign(sample=s) for s, df in feat_data.items()], ignore_index=True)
    feat_sub = (feat_plot_df
                .groupby(['sample', 'feature_type'], group_keys=False)
                .apply(lambda g: g.sample(min(len(g), 50_000), random_state=42))
                .reset_index(drop=True))

    FT_ORDER   = ['promoter', 'gene_body', 'intergenic']
    FT_LABELS  = {'promoter': f'Promoter\n({promoter_bp // 1000} kb)',
                  'gene_body': 'Gene body', 'intergenic': 'Intergenic'}
    FT_PALETTE = {'promoter': '#9C27B0', 'gene_body': '#FF9800',
                  'intergenic': '#607D8B'}
    sorder = [s for s in sample_order if s in feat_data]

    fig, axes = plt.subplots(1, len(sorder), figsize=(4 * len(sorder), 5),
                             sharey=True)
    if len(sorder) == 1:
        axes = [axes]
    for ax, sample in zip(axes, sorder):
        sub = feat_sub[feat_sub['sample'] == sample]
        sns.violinplot(data=sub, x='feature_type', y='pct_mod', order=FT_ORDER,
                       palette=FT_PALETTE, inner='quartile', linewidth=0.8, ax=ax)
        ax.set_title(SAMPLE_META.get(sample, {}).get('label', sample))
        ax.set_xlabel('')
        ax.set_ylabel('CpG methylation (%)' if ax == axes[0] else '')
        ax.set_xticklabels([FT_LABELS[ft] for ft in FT_ORDER], fontsize=8)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter())
        ax.set_ylim(-5, 105)
        sns.despine(ax=ax)
    fig.suptitle('CpG methylation by genomic feature type', y=1.02, fontsize=13)
    plt.tight_layout()
    out = out_dir / 'fig2_methylation_by_feature.pdf'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log.info(f'  Saved: {out}')

    return feat_data


# ─────────────────────────────────────────────────────────────────────────────
# Section 3 — PCA and hierarchical clustering
# ─────────────────────────────────────────────────────────────────────────────

def section3(results_dir, available, out_dir, min_cov=5, top_n_sites=50_000):
    log.info('=== Section 3: PCA / hierarchical clustering ===')
    cpg_samples = available.get('CpG', [])
    if len(cpg_samples) < 2:
        log.warning('Fewer than 2 CpG samples available — skipping Section 3.')
        return

    frames = []
    for sample in cpg_samples:
        fpath = results_dir / 'CpG' / sample / f'{sample}.bedMethyl'
        log.info(f'  Loading {sample} ...')
        bm = load_bedmethyl(fpath, usecols=[0, 1, 9, 10])
        bm = bm[['chrom', 'start', 'N_valid_cov', 'pct_mod']]
        bm = bm[bm['N_valid_cov'] >= min_cov]
        bm = bm.rename(columns={'pct_mod': sample})[['chrom', 'start', sample]]
        frames.append(bm)
        log.info(f'    {len(bm):,} sites (cov ≥ {min_cov})')

    matrix = reduce(
        lambda a, b: a.merge(b, on=['chrom', 'start'], how='inner'), frames)
    log.info(f'  Common sites (all samples, cov ≥ {min_cov}): {len(matrix):,}')

    site_var  = matrix[cpg_samples].var(axis=1)
    top_sites = site_var.nlargest(min(top_n_sites, len(matrix))).index
    mat_filt  = matrix.loc[top_sites, cpg_samples]
    log.info(f'  Top {len(mat_filt):,} most variable sites retained')

    X = mat_filt.T.values     # (n_samples × n_sites)

    # PCA
    pca   = PCA(n_components=min(len(cpg_samples), 3))
    X_sc  = StandardScaler().fit_transform(X)
    X_pca = pca.fit_transform(X_sc)
    var_exp = pca.explained_variance_ratio_ * 100
    log.info('  Explained variance: ' +
             '  '.join(f'PC{i+1}: {v:.1f}%' for i, v in enumerate(var_exp)))

    fig, ax = plt.subplots(figsize=(6, 5))
    for i, sample in enumerate(cpg_samples):
        ax.scatter(X_pca[i, 0], X_pca[i, 1],
                   color=SAMPLE_META.get(sample, {}).get('color', 'grey'),
                   s=160, zorder=3, edgecolors='white', linewidth=0.8)
        ax.annotate(SAMPLE_META.get(sample, {}).get('label', sample),
                    (X_pca[i, 0], X_pca[i, 1]),
                    textcoords='offset points', xytext=(8, 4), fontsize=9)
    ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}%)' if pca.n_components_ > 1 else 'PC2')
    ax.set_title(f'PCA of CpG methylation\n'
                 f'(top {len(mat_filt):,} variable sites, cov ≥ {min_cov})')
    ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')
    sns.despine(ax=ax)
    plt.tight_layout()
    out = out_dir / 'fig3a_pca.pdf'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log.info(f'  Saved: {out}')

    # Clustermap
    rng    = np.random.default_rng(0)
    hm_idx = rng.choice(len(mat_filt),
                        min(10_000, len(mat_filt)), replace=False)
    hm_mat = mat_filt.iloc[hm_idx]

    cg = sns.clustermap(
        hm_mat.T,
        col_cluster=True,
        row_cluster=(len(cpg_samples) > 2),
        cmap='RdBu_r', center=50, vmin=0, vmax=100,
        xticklabels=False,
        yticklabels=[SAMPLE_META.get(s, {}).get('label', s)
                     for s in hm_mat.columns],
        figsize=(12, max(2 * len(cpg_samples), 5)),
        cbar_kws={'label': 'CpG methylation (%)'},
    )
    cg.fig.suptitle(
        f'Hierarchical clustering — CpG methylation\n'
        f'({len(hm_idx):,} random sites from top {len(mat_filt):,} variable)',
        y=1.02,
    )
    out = out_dir / 'fig3b_clustermap.pdf'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log.info(f'  Saved: {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Section 4 — Coverage and methylation distribution (ONT QC)
# ─────────────────────────────────────────────────────────────────────────────

def section4(allC_data, allC_ctx, out_dir, sample_order):
    log.info('=== Section 4: Coverage and methylation distribution ===')
    if not allC_data:
        log.warning('No allC data loaded — skipping Section 4.')
        return

    sorder = [s for s in sample_order if s in allC_data]

    # ── fig4a: coverage histogram + per-chromosome box ────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    ax = axes[0]
    for s in sorder:
        cov = allC_data[s]['N_valid_cov'].clip(upper=200)
        ax.hist(cov, bins=80, alpha=0.5,
                label=SAMPLE_META.get(s, {}).get('label', s),
                color=SAMPLE_META.get(s, {}).get('color', None), edgecolor='none')
    ax.set_xlabel('Coverage (N_valid_cov, clipped at 200×)')
    ax.set_ylabel('Number of cytosine positions')
    ax.set_title('Coverage distribution (allC)')
    ax.legend(frameon=False, fontsize=8)
    sns.despine(ax=ax)

    ax = axes[1]
    chrom_frames = []
    for s in sorder:
        by_chr = (allC_data[s].groupby('chrom')['N_valid_cov']
                              .median().reset_index())
        by_chr['sample'] = s
        chrom_frames.append(by_chr)
    chrom_cov_df = pd.concat(chrom_frames, ignore_index=True)
    chrom_cov_df = chrom_cov_df[chrom_cov_df['chrom'].str.startswith('NC_')]
    sns.boxplot(data=chrom_cov_df, x='sample', y='N_valid_cov', order=sorder,
                palette=[SAMPLE_META.get(s, {}).get('color', 'grey') for s in sorder],
                width=0.5, linewidth=0.8, ax=ax)
    ax.set_xticklabels(
        [SAMPLE_META.get(s, {}).get('label', s) for s in sorder])
    ax.set_xlabel('')
    ax.set_ylabel('Median coverage per chromosome')
    ax.set_title('Per-chromosome median coverage (main chromosomes)')
    sns.despine(ax=ax)
    plt.tight_layout()
    out = out_dir / 'fig4a_coverage_distribution.pdf'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log.info(f'  Saved: {out}')

    # ── fig4b: methylation % distribution per context × sample ────────────────
    ctx_order = ['CG', 'CHG', 'CHH']
    n_ctx     = len(ctx_order)
    n_sam     = len(sorder)
    fig, axes = plt.subplots(n_ctx, n_sam, figsize=(4 * n_sam, 3 * n_ctx),
                             sharex=True, sharey='row')
    if n_sam == 1:
        axes = axes.reshape(-1, 1)
    for ri, ctx in enumerate(ctx_order):
        for ci, s in enumerate(sorder):
            ax = axes[ri][ci]
            if s not in allC_ctx:
                ax.set_visible(False)
                continue
            sub = allC_ctx[s][allC_ctx[s]['context'] == ctx]
            ax.hist(sub['pct_mod'], bins=50, edgecolor='none', alpha=0.8,
                    color=SAMPLE_META.get(s, {}).get('color', None))
            ax.set_title(
                f'{SAMPLE_META.get(s, {}).get("label", s)} — {ctx}', fontsize=9)
            if ci == 0:
                ax.set_ylabel('Count')
            if ri == n_ctx - 1:
                ax.set_xlabel('Methylation (%)')
            ax.xaxis.set_major_formatter(mticker.PercentFormatter())
            sns.despine(ax=ax)
    fig.suptitle('Methylation % distribution by context and sample (allC)', y=1.01)
    plt.tight_layout()
    out = out_dir / 'fig4b_methylation_distribution.pdf'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log.info(f'  Saved: {out}')

    # ── fig4c: hexbin pct_mod vs coverage ─────────────────────────────────────
    fig, axes = plt.subplots(1, n_sam, figsize=(5 * n_sam, 4), sharey=True)
    if n_sam == 1:
        axes = [axes]
    for ax, s in zip(axes, sorder):
        df = allC_data[s]
        hb = ax.hexbin(df['N_valid_cov'].clip(upper=150), df['pct_mod'],
                       gridsize=50, cmap='YlOrRd', mincnt=1, bins='log')
        ax.set_xlabel('Coverage (clipped at 150×)')
        ax.set_ylabel('Methylation (%)' if ax == axes[0] else '')
        ax.set_title(SAMPLE_META.get(s, {}).get('label', s))
        ax.yaxis.set_major_formatter(mticker.PercentFormatter())
        plt.colorbar(hb, ax=ax, label='log10(count)')
        sns.despine(ax=ax)
    fig.suptitle('Methylation % vs coverage — hexbin (allC, all contexts)',
                 y=1.02, fontsize=12)
    plt.tight_layout()
    out = out_dir / 'fig4c_methylation_vs_coverage_hexbin.pdf'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log.info(f'  Saved: {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Methylation landscape EDA (headless, all 4 sections).')
    p.add_argument('--wdir', default='/home/daffa/Work/2026/02-JSPP67',
                   help='Project root directory (default: %(default)s)')
    p.add_argument('--ref-fasta', default=None,
                   help='Path to reference FASTA (default: wdir/data/reference/asm_BTx623.fna)')
    p.add_argument('--ref-gff', default=None,
                   help='Path to annotation GFF3 (default: wdir/data/reference/GCF_*_genomic.gff)')
    p.add_argument('--results-dir', default=None,
                   help='DMR results directory (default: wdir/results/dmr_analysis)')
    p.add_argument('--out-dir', default=None,
                   help='Output directory for figures (default: wdir/results/methylation_landscape)')
    p.add_argument('--samples', nargs='+',
                   default=['SBC4', 'SBC10', 'SBC11', 'SBC23'],
                   help='Sample names (default: SBC4 SBC10 SBC11 SBC23)')
    p.add_argument('--contexts', nargs='+', default=['CpG', 'allC'],
                   help='Methylation contexts to consider (default: CpG allC)')
    p.add_argument('--min-cov', type=int, default=5,
                   help='Minimum coverage for PCA matrix (default: 5)')
    p.add_argument('--top-n-sites', type=int, default=50_000,
                   help='Most variable CpG sites to retain for PCA (default: 50000)')
    p.add_argument('--promoter-bp', type=int, default=2000,
                   help='Promoter window upstream of TSS in bp (default: 2000)')
    return p.parse_args()


def main():
    args = parse_args()

    wdir        = Path(args.wdir)
    ref_fasta   = Path(args.ref_fasta)  if args.ref_fasta   else wdir / 'data/reference/asm_BTx623.fna'
    ref_gff     = Path(args.ref_gff)    if args.ref_gff     else wdir / 'data/reference/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff'
    results_dir = Path(args.results_dir) if args.results_dir else wdir / 'results/dmr_analysis'
    out_dir     = Path(args.out_dir)    if args.out_dir     else wdir / 'results/methylation_landscape'
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info(f'WDIR        : {wdir}')
    log.info(f'REF_FASTA   : {ref_fasta}  exists={ref_fasta.exists()}')
    log.info(f'REF_GFF     : {ref_gff}  exists={ref_gff.exists()}')
    log.info(f'RESULTS_DIR : {results_dir}')
    log.info(f'OUT_DIR     : {out_dir}')
    log.info(f'SAMPLES     : {args.samples}')
    log.info(f'CONTEXTS    : {args.contexts}')

    available = discover_available(results_dir, args.samples, args.contexts)
    for ctx, slist in available.items():
        log.info(f'  {ctx}: {slist}')

    tmpdir = Path(tempfile.mkdtemp(prefix='methyl_eda_'))
    log.info(f'Temporary directory: {tmpdir}')

    # Section 1
    s1_result = section1(results_dir, available, ref_fasta, out_dir, args.samples)
    allC_data = s1_result[0] if isinstance(s1_result, tuple) else {}
    allC_ctx  = s1_result[1] if isinstance(s1_result, tuple) else {}

    # Section 2
    section2(results_dir, available, ref_gff, out_dir, tmpdir,
             args.samples, promoter_bp=args.promoter_bp)

    # Section 3
    section3(results_dir, available, out_dir,
             min_cov=args.min_cov, top_n_sites=args.top_n_sites)

    # Section 4
    section4(allC_data, allC_ctx, out_dir, args.samples)

    log.info('Done. All figures written to: %s', out_dir)


if __name__ == '__main__':
    main()
