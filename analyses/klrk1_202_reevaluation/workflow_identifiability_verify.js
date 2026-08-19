export const meta = {
  name: 'klrk1-202-203-identifiability-verify',
  description: 'Adversarially verify the exon-overlap / EM-identifiability claims underlying the Klrk1-202 exclusion decision and the trust we place in Klrk1-203',
  phases: [
    { title: 'Verify', detail: 'one adversarial verifier per structural/read-level claim' },
    { title: 'Synthesize', detail: 'calibrated trust statement + list of my claims needing correction' },
  ],
}

const CTX = `
CONTEXT — repo: /Users/study/Desktop/Karimi/Klrk1_GVHD_project  (cd here first)
Tools available: bash, python3 (numpy/pandas), samtools at /opt/homebrew/bin/samtools.
GTF with Klrk1 exons: shared/ref/genome/Klrk1.gtf  (exon lines: col3==exon, col4=start,col5=end,col7=strand, attrs col9 has transcript_id "ENSMUST....").
Transcript IDs -> names: 201=ENSMUST00000032252, 202=ENSMUST00000095412, 203=ENSMUST00000137660,
204=ENSMUST00000152256, 205=ENSMUST00000168919, 206=ENSMUST00000204694.
Biotype: 201/202/205 protein_coding; 203/204/206 retained_intron.
Per-sample 6-isoform TPM table: python_pipeline/outputs/tables/Supplementary_all_sample_6isoform.csv
(cols: dataset,condition,sample,TPM_201,TPM_202,TPM_203,TPM_204,TPM_205,TPM_206,total6_TPM).
BAMs (chrom token = "6", not "chr6"): datasets/GSE147371_GVHD_CD4/data/bam/GVHD_CD4_Tem_1.bam, _Tem_2.bam, _Tn.bam ;
datasets/GSE203167_Karimi_TCF7/data/bam/WT_Pre.bam, WT_Post7.bam. Index with samtools index if no .bai.
Key coordinates (GRCm39 chr6, minus strand):
  - 203 intron-4-proximal retention segment: 129593531-129593631 (101bp)
  - 204-unique retained region (intron-3, 204 retains / 203 splices): 129593737-129594445 (709bp)
  - constitutive exon4: 129593632-129593736
  - 202/205 alternative first exon start: 129600774 (202 ends 129600804; 205 ends 129600827)
  - 201 canonical first exon: 129599500-129599735
You MUST run commands to check; do not answer from priors. Reproduce numbers yourself.
Return your verdict as: CONFIRMED (claim fully holds), PARTIAL (holds with an important caveat/correction),
or REFUTED (claim is wrong). Always give the exact numbers you obtained.
`

const SCHEMA = {
  type: 'object',
  properties: {
    claim_id: { type: 'string' },
    verdict: { type: 'string', enum: ['CONFIRMED', 'PARTIAL', 'REFUTED'] },
    numbers_obtained: { type: 'string', description: 'the concrete values you computed' },
    reasoning: { type: 'string' },
    correction: { type: 'string', description: 'if PARTIAL/REFUTED, the corrected statement; else empty' },
  },
  required: ['claim_id', 'verdict', 'numbers_obtained', 'reasoning', 'correction'],
}

const CLAIMS = [
  { id: 'C1_unique_bp', prompt: `CLAIM C1: Among the six Klrk1 transcripts, 201, 202, and 203 each have ZERO transcript-exclusive exonic bases (every base is shared with >=1 other Klrk1 transcript), whereas 204 (~784bp), 205 (~2113bp), and 206 (~2865bp) have large exclusive regions. Independently recompute per-transcript exonic length and exclusive-bp (bases in this transcript's exons not in the union of all other Klrk1 transcripts' exons) by parsing the GTF. Try to REFUTE.` },
  { id: 'C2_202_subset_205', prompt: `CLAIM C2: Klrk1-202 is a strict subsequence of Klrk1-205 (every exonic base of 202 is also exonic in 205), so NO read can distinguish 202 from 205 and the 202-vs-205 Salmon split is fundamentally unidentifiable. Independently compute |bases(202) \\ bases(205)| from the GTF. Also compute |bases(202)\\bases(201)| (should be ~31, the alt first exon). Try to REFUTE.` },
  { id: 'C3_203_vs_204_structure', prompt: `CLAIM C3: 203's intron-4-proximal retention segment (129593531-129593631) is shared ONLY with 204 (not 201/202/205/206), AND 203 is read-level distinguishable from 204 because 204 additionally retains intron-3 (129593737-129594445, 709bp, 204-exclusive) which 203 splices out. Verify from the GTF: (a) which transcripts cover 129593531-129593631; (b) that 129593737-129594445 is 204-exclusive. Try to REFUTE.` },
  { id: 'C4_bam_203_induced', prompt: `CLAIM C4: In GVHD Tem BAMs, mean depth over the 203+204 retention segment (129593531-129593631) is ~20-24x and EXCEEDS the depth over the 204-exclusive intron-3 region (129593737-129594445, ~10x); in naive Tn both are ~0 (retention 0.4, exon4 9). Therefore a 203-type contribution (retains intron-4-proximal but splices intron-3) is required and is differentiation-induced, not merely 204. Recompute samtools depth for Tem_1, Tem_2, Tn over: 129593531-129593631, 129593737-129594445, 129593632-129593736 (exon4). ADVERSARIAL: is the excess retention over 204's level genuinely attributable to 203, or could it be explained otherwise (e.g. 206, coverage bias, or the segment also being intronic in 204)? State the strongest counter-argument and whether it survives.` },
  { id: 'C5_correlation', prompt: `CLAIM C5: 201 and 202 POSITIVELY correlate across samples (NOT anticorrelated); within GSE203167 r(201,202)~+0.97; and pooling 201+202 does NOT consistently reduce coefficient-of-variation vs 202 alone. (This corrects an earlier claim that pooling reduces CV via anticorrelation.) Using Supplementary_all_sample_6isoform.csv, compute Pearson r(201,202), r(205,202) across samples where TPM_202>0, and within GSE203167; and compare within-condition CV of 202 vs (201+202) vs (202+205) for GSE203167 conditions and GSE119943 Arm_D4.5/pMIG. Report where pooling helps and where it does not. Try to REFUTE the corrected claim.` },
  { id: 'C6_denominator', prompt: `CLAIM C6: Because the reads Salmon assigns to 202 are real Klrk1 reads (they map to Klrk1 exons; 202 shares the alt first exon 129600774 with 205/206 which have genuine BAM coverage ~10-13x in GSE203167 WT), 202 must be counted in the total-Klrk1 denominator; excluding only 202 while keeping the near-identical 205 undercounts total Klrk1 and inflates RI% and 203%. Verify the alt-first-exon coverage in GSE203167 WT_Pre/WT_Post7 BAMs (129600774-129600804), and confirm 202 being protein_coding means it enters the denominator but not the RI numerator. Try to REFUTE the logic that 202 belongs in the denominator.` },
]

phase('Verify')
const verdicts = await parallel(CLAIMS.map(c => () =>
  agent(CTX + '\n\n' + c.prompt, { label: c.id, phase: 'Verify', schema: SCHEMA, effort: 'high' })
    .then(v => ({ ...v, claim_id: v?.claim_id || c.id }))
))
const clean = verdicts.filter(Boolean)

phase('Synthesize')
const synth = await agent(
  CTX + '\n\nYou are given adversarial verification verdicts (JSON) on 6 claims about Klrk1 isoform ' +
  'identifiability:\n' + JSON.stringify(clean, null, 2) +
  '\n\nWrite a calibrated synthesis answering the user\'s core questions: ' +
  '(1) How much should we TRUST the Klrk1-203 fraction, given 203 has 0 exclusive bp and shares its ' +
  'retention with 204? Distinguish "203 induction DIRECTION" (robustness) from "exact 203% vs 204" (model-dependence). ' +
  '(2) Why is 202 in a WORSE position than 203 (202 subset of 205 vs 203 partially identifiable from 204)? ' +
  '(3) Which of the 6 claims were CONFIRMED / PARTIAL / REFUTED, and list any statement in my prior ' +
  'decision note that must be corrected (esp. any anticorrelation/pooling-CV claim). ' +
  'Be adversarial and precise; if the evidence does not support a claim, say so.',
  { label: 'synthesize', phase: 'Synthesize', effort: 'high' })

return { verdicts: clean, synthesis: synth }
