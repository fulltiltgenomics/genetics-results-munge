# rCNV sliding windows: GRCh37 -> GRCh38 liftOver measurement

Collins et al. 2022 (Cell 185:3041, Zenodo 6347673) ships `sliding_window_sumstats.tar.gz`:
108 bgzipped BED files (54 HPO phenotypes x DEL/DUP) over **267,237 autosomal 200 kb windows
in 10 kb steps, GRCh37**. The suite is GRCh38 throughout, so the windows either lift or the
product is dropped — the epic's second descope trigger is *more than ~2% of windows failing to
lift cleanly*.

This document is that measurement. **Every number below is pasted from
`scripts/rcnv_liftover_windows.py`'s stdout, not typed** — re-run the script to reproduce it.

The measurement is not just a record: `scripts/munge_rcnv.py --product windows` performs the
same lift over the same window set and asserts it reproduces the lifted/dropped split below
exactly, so a chain or binary that lifts differently stops the munge rather than quietly
shipping a different window set. See `docs/rcnv-dosage-sensitivity.md` for what that product
then does with the lifted windows.

## DECISION

**1.826% of the 267,237 windows fail to lift cleanly, below the ~2% trigger, so the sliding-window
product SHIPS and hpa1.12 proceeds.** 262,357 windows (98.174%) map to the same chromosome with a
lifted length inside 180-220 kb.

The margin is thin (1.826% against 2.000%), and it is thin for a reason that is not noise: 4,880 of
the failures are concentrated in a few dozen GRCh37 regions — chr9 alone loses 7.685% of its windows,
almost all to the pericentromeric/heterochromatic block, and chr21 (4.358%), chr22 (3.616%) and
chr1 (3.242%) each lose several percent (per-chromosome table below).

The shipped tolerance (lifted length in 180-220 kb, UCSC-default `minMatch` 0.95) was pre-registered
by the task before this measurement ran, not chosen after seeing the number. Measured flip points:
widening to ±5% (190-210 kb) gives 1.853% (still ships); narrowing to ±2.5% (195-205 kb) gives
2.007% (dropped). `minMatch=0.99` gives 1.948% unmapped before the chromosome/length filters and
would exceed 2% after them; `minMatch=1.0` gives 7.048%. The margin between the shipped 1.826% and
the 2% trigger is 465 windows. The tolerance is recorded in the header of every run and should not
be changed without re-recording this decision.

## Command

```
python3 scripts/rcnv_liftover_windows.py \
    --sumstats-dir <unpacked>/Collins_rCNV_2022.sliding_window_sumstats \
    --liftover-bin <dir>/liftOver \
    --chain <dir>/hg19ToHg38.over.chain.gz \
    --work-dir <dir>/work
```

Tool and chain, both fetched from UCSC for this run (neither is vendored in the repo):

| what | source | identity |
|---|---|---|
| chain | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz` | md5 `35887f73fe5e2231656504d1f6430900`, `Last-Modified: 2014-01-01` |
| binary | `https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver` | `Last-Modified: 2026-08-20`, 30,663,104 bytes |

## Procedure

Taken from `scripts/munge_calderon.py`, this repo's only liftOver munge (`liftover_intervals()`):
the **whole interval** is lifted as one BED4 record whose name is the GRCh37 key, with UCSC
defaults (`minMatch` 0.95, no `-multiple`); keys that map more than once are dropped as
multi-mapped, and keys landing on a different chromosome than their GRCh37 chromosome are dropped
as inconsistent. Calderon's peaks have no expected width so it stops there; a 200 kb window does,
so this measurement adds the length filter — lifted length must be within **180-220 kb**.

liftOver reports an interval as unmapped without saying which end failed, so a second pass lifts
each window's two 1 bp endpoints independently. That pass only *attributes* an interval failure to
an end; it never rescues a window. A window whose two endpoints both lift but whose interval
liftOver still refused is `interior_deleted` — liftOver's "Partially deleted in new", i.e. more
than 5% of the span has no hg38 image.

Failure classes, in the order the report prints them:

| class | meaning |
|---|---|
| `unmapped_start` | interval unmapped; the 5' endpoint has no hg38 image |
| `unmapped_end` | interval unmapped; the 3' endpoint has no hg38 image |
| `unmapped_both` | interval unmapped; neither endpoint has an hg38 image |
| `split` | liftOver "Split in new" (or the key mapped more than once) |
| `chromosome_changed` | lifted to a different chromosome than the GRCh37 one |
| `length_out_of_tolerance` | same chromosome, lifted length outside 180-220 kb |
| `interior_deleted` | both endpoints lift, but >5% of the interval has no hg38 image |

## What hpa1.12 has to carry

- **Lifted length is not constant.** 90.541% of surviving windows are exactly 200,000 bp; the rest
  range 190,000-219,265. A consumer must not assume a 200 kb width from the lifted coordinates —
  which is why the table stores `start_grch37`/`end_grch37` alongside.
- **Neighbour geometry mostly survives, but not everywhere.** Of 262,335 adjacent same-chromosome
  pairs, 378 (0.144%) come out *reordered* — the later GRCh37 window starts at or before its
  predecessor in GRCh38 — and a further 218 (0.083%) overlap by more than one 10 kb step away from
  their GRCh37 overlap. 758 surviving windows (0.289%) sit in such a pair. A "next window" or
  "adjacent window" query answered by GRCh38 ordering will disagree with GRCh37 ordering for those.
- **The dropped windows are not spread evenly**, so a region query over a lost block returns
  nothing rather than fewer rows. The per-chromosome and per-region tables below are the record of
  which blocks those are; the full per-window list is written to `dropped_windows.tsv` in the run's
  `--work-dir`.

## Run output (verbatim)

```
==============================================================================
rCNV sliding windows — UCSC liftOver hg19 -> hg38
==============================================================================
liftOver binary:           /tmp/claude-1001/-home-jkarjala-staging-genetics-results-suite/a9dd2f16-2506-4ede-83e9-573d341b2e8c/scratchpad/liftover/liftOver
chain:                     hg19ToHg38.over.chain.gz  md5 35887f73fe5e2231656504d1f6430900
chain source:              https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
tolerance:                 lifted length in [180,000, 220,000] bp, same chromosome

sumstats files:            108
window set verified on:    2 files (HP0000118.DEL, HP0000118.DUP)
  sha256(chr,start,end)    1bada03bc8ab01de  HP0000118.rCNV.DEL.sliding_window.meta_analysis.stats.bed.gz
  sha256(chr,start,end)    1bada03bc8ab01de  HP0000118.rCNV.DUP.sliding_window.meta_analysis.stats.bed.gz
window set identical:      yes (267237 windows per file)

liftOver reasons reported for intervals not mapped (raw, pass 1):
  Split in new                          2,382
  Partially deleted in new              2,051

WINDOW COUNTS
  total windows                       267,237  100.000%
  mapped cleanly                      262,357   98.174%
    failed: unmapped_start                576    0.216%
    failed: unmapped_end                  608    0.228%
    failed: unmapped_both                  34    0.013%
    failed: split                       2,382    0.891%
    failed: chromosome_changed              4    0.001%
    failed: length_out_of_tolerance       443    0.166%
    failed: interior_deleted              833    0.312%
  failed (all classes)                  4,880    1.826%

LIFTED LENGTH (clean windows, bp)
  min                                 190,000
  median                              200,000
  max                                 219,265
  exactly 200,000                     237,541   90.541%

NEIGHBOUR GEOMETRY (adjacent clean windows on the same chromosome, 262,335 pairs)
  pairs reordered after lifting           378    0.144%
  pairs whose overlap moved > step        218    0.083%
  clean windows in such a pair            758    0.289%

DROPPED WINDOWS BY CHROMOSOME
  chr   windows  dropped       %    unm_s   unm_e  unm_se   split chr_chg len_oot int_del
    1    22,514      730  3.242%       62      57       6     515       4      40      46
    2    23,681      269  1.136%       31      40       0      91       0       9      98
    3    19,456       71  0.365%       15      15       0      21       0      20       0
    4    18,760      159  0.848%       23      29       2      23       0       0      82
    5    17,764      103  0.580%       21      21       0       1       0       0      60
    6    16,227      227  1.399%       26      36       0      85       0      65      15
    7    15,550      367  2.360%       30      48       2     174       0      40      73
    8    14,268      159  1.114%       18      23       0      83       0       5      30
    9    11,997      922  7.685%       70      73      15     687       0      17      60
   10    13,128      403  3.070%       45      34       0     258       0       6      60
   11    13,100      171  1.305%       30      23       0      47       0      40      31
   12    13,039      104  0.798%       17      22       0      20       0       0      45
   13     9,518       51  0.536%       26      25       0       0       0       0       0
   14     8,691       69  0.794%        4       0       0      25       0      40       0
   15     8,104      180  2.221%       13      13       0     125       0       6      23
   16     7,837      114  1.455%       18      23       0       8       0      35      30
   17     7,756      198  2.553%       23      22       0      63       0      60      30
   18     7,447       91  1.222%       18      18       0       5       0      20      30
   19     5,574      110  1.973%       13      13       0      80       0       0       4
   20     5,941      107  1.801%       27      28       0       0       0      20      32
   21     3,511      153  4.358%       32      28       9      15       0       0      69
   22     3,374      122  3.616%       14      17       0      56       0      20      15

DROPPED WINDOWS BY REGION (contiguous runs of one class; full list: dropped_windows.tsv)
  chr  GRCh37 region              class                    windows
    1  0-460,000                  split                         27
    1  270,000-470,000            unmapped_start                 1
    1  280,000-510,000            unmapped_both                  4
    1  320,000-520,000            unmapped_end                   1
    1  330,000-670,000            split                         15
    1  480,000-710,000            unmapped_start                 4
    1  2,450,000-2,680,000        unmapped_end                   4
    1  2,490,000-2,830,000        interior_deleted              15
    1  2,640,000-2,870,000        unmapped_start                 4
    1  3,660,000-3,900,000        unmapped_end                   5
    1  3,940,000-4,180,000        unmapped_start                 5
    1  12,830,000-13,890,000      split                         69
    1  16,940,000-17,170,000      unmapped_end                   4
    1  16,980,000-17,320,000      interior_deleted              15
    1  17,130,000-17,360,000      unmapped_start                 4
    1  29,690,000-29,930,000      unmapped_end                   5
    1  29,970,000-30,210,000      unmapped_start                 5
    1  103,600,000-104,060,000    split                         27
    1  103,870,000-104,100,000    unmapped_start                 4
    1  120,510,000-120,990,000    split                         29
    1  121,030,000-121,270,000    unmapped_start                 5
    1  121,300,000-121,540,000    unmapped_end                   5
    1  142,480,000-142,720,000    unmapped_start                 5
    1  142,530,000-142,740,000    chromosome_changed             2
    1  142,550,000-142,780,000    unmapped_end                   4
    1  142,590,000-142,930,000    split                         15
    1  142,740,000-142,960,000    unmapped_start                 3
    1  142,770,000-142,980,000    unmapped_both                  2
    1  142,790,000-143,020,000    unmapped_end                   4
    1  143,060,000-143,490,000    split                         24
    1  143,300,000-143,530,000    unmapped_start                 4
    1  143,340,000-143,550,000    chromosome_changed             2
    1  143,360,000-143,600,000    unmapped_end                   5
    1  143,590,000-145,000,000    split                         86
    1  144,980,000-145,370,000    length_out_of_tolerance       20
    1  145,650,000-145,880,000    unmapped_end                   4
    1  145,690,000-146,640,000    split                         76
    1  147,750,000-148,080,000    split                         14
    1  148,120,000-148,420,000    split                         11
    1  148,460,000-149,220,000    split                         57
    1  149,270,000-149,480,000    split                          2
    1  149,290,000-149,500,000    unmapped_end                   2
    1  149,310,000-149,510,000    interior_deleted               1
    1  149,320,000-149,820,000    split                         31
    1  205,740,000-205,980,000    unmapped_end                   5
    1  206,020,000-206,260,000    unmapped_start                 5
    1  206,150,000-206,390,000    split                          5
    1  206,430,000-206,670,000    unmapped_start                 5
    1  223,540,000-223,940,000    split                         21
    1  223,750,000-223,980,000    unmapped_start                 4
    1  228,550,000-228,940,000    length_out_of_tolerance       20
    1  235,010,000-235,240,000    unmapped_end                   4
    1  235,050,000-235,390,000    interior_deleted              15
    1  235,200,000-235,430,000    unmapped_start                 4
    1  248,720,000-248,960,000    unmapped_end                   5
    1  249,000,000-249,250,000    split                          6
    2  3,340,000-3,570,000        unmapped_end                   4
    2  3,380,000-3,720,000        interior_deleted              15
    2  3,530,000-3,760,000        unmapped_start                 4
    2  4,830,000-5,070,000        unmapped_end                   5
    2  5,060,000-5,300,000        unmapped_start                 5
    2  16,090,000-16,320,000      unmapped_end                   4
    2  16,130,000-16,470,000      interior_deleted              15
    2  16,280,000-16,510,000      unmapped_start                 4
    2  19,020,000-19,410,000      split                         20
    2  20,970,000-21,170,000      unmapped_end                   1
    2  20,980,000-21,350,000      interior_deleted              18
    2  21,160,000-21,360,000      unmapped_start                 1
    2  87,470,000-87,910,000      split                         25
    2  91,540,000-91,910,000      split                         18
    2  92,140,000-92,380,000      unmapped_end                   5
    2  95,270,000-95,510,000      split                          5
    2  97,830,000-98,210,000      split                         19
    2  98,020,000-98,300,000      length_out_of_tolerance        9
    2  109,920,000-110,160,000    unmapped_end                   5
    2  110,200,000-110,440,000    unmapped_start                 5
    2  149,510,000-149,750,000    unmapped_end                   5
    2  149,740,000-149,980,000    unmapped_start                 5
    2  233,820,000-234,050,000    unmapped_end                   4
    2  233,860,000-234,200,000    interior_deleted              15
    2  234,010,000-234,240,000    unmapped_start                 4
    2  239,620,000-239,830,000    unmapped_end                   2
    2  239,640,000-240,000,000    interior_deleted              17
    2  239,810,000-240,020,000    unmapped_start                 2
    2  240,600,000-240,800,000    unmapped_end                   1
    2  240,610,000-240,980,000    interior_deleted              18
    2  240,790,000-240,990,000    unmapped_start                 1
    2  242,920,000-243,150,000    unmapped_end                   4
    2  242,960,000-243,190,000    split                          4
    3  10,000-250,000             unmapped_start                 5
    3  57,190,000-57,580,000      length_out_of_tolerance       20
    3  65,990,000-66,230,000      unmapped_end                   5
    3  66,220,000-66,460,000      unmapped_start                 5
    3  90,320,000-90,560,000      unmapped_end                   5
    3  93,450,000-93,690,000      unmapped_start                 5
    3  195,010,000-195,410,000    split                         21
    3  197,780,000-198,020,000    unmapped_end                   5
    4  0-200,000                  unmapped_start                 1
    4  1,240,000-1,470,000        unmapped_end                   4
    4  1,280,000-1,620,000        interior_deleted              15
    4  1,430,000-1,660,000        unmapped_start                 4
    4  8,610,000-8,810,000        unmapped_end                   1
    4  8,620,000-8,990,000        interior_deleted              18
    4  8,800,000-9,000,000        unmapped_start                 1
    4  9,090,000-9,320,000        unmapped_end                   4
    4  9,130,000-9,470,000        interior_deleted              15
    4  9,280,000-9,510,000        unmapped_start                 4
    4  31,640,000-32,020,000      interior_deleted              19
    4  49,150,000-49,390,000      unmapped_end                   5
    4  49,430,000-49,660,000      unmapped_start                 4
    4  49,470,000-49,680,000      unmapped_both                  2
    4  49,490,000-49,720,000      unmapped_end                   4
    4  52,610,000-52,850,000      unmapped_start                 5
    4  59,550,000-59,780,000      unmapped_end                   4
    4  59,590,000-59,930,000      interior_deleted              15
    4  59,740,000-59,970,000      unmapped_start                 4
    4  75,240,000-75,450,000      unmapped_end                   2
    4  75,260,000-75,680,000      split                         23
    4  190,860,000-191,100,000    unmapped_end                   5
    5  0-200,000                  split                          1
    5  17,350,000-17,580,000      unmapped_end                   4
    5  17,390,000-17,730,000      interior_deleted              15
    5  17,540,000-17,770,000      unmapped_start                 4
    5  46,220,000-46,460,000      unmapped_end                   5
    5  49,350,000-49,590,000      unmapped_start                 5
    5  91,450,000-91,680,000      unmapped_end                   4
    5  91,490,000-91,830,000      interior_deleted              15
    5  91,640,000-91,870,000      unmapped_start                 4
    5  138,600,000-138,830,000    unmapped_end                   4
    5  138,640,000-138,980,000    interior_deleted              15
    5  138,790,000-139,020,000    unmapped_start                 4
    5  154,950,000-155,180,000    unmapped_end                   4
    5  154,990,000-155,330,000    interior_deleted              15
    5  155,140,000-155,370,000    unmapped_start                 4
    6  10,000-240,000             unmapped_start                 4
    6  26,550,000-26,950,000      split                         21
    6  26,760,000-26,990,000      length_out_of_tolerance        4
    6  50,900,000-51,290,000      length_out_of_tolerance       20
    6  57,120,000-57,510,000      length_out_of_tolerance       20
    6  57,900,000-58,130,000      unmapped_end                   4
    6  57,940,000-58,280,000      split                         15
    6  58,090,000-58,320,000      unmapped_start                 4
    6  58,140,000-58,340,000      length_out_of_tolerance        1
    6  58,600,000-58,840,000      unmapped_end                   5
    6  61,830,000-62,070,000      split                          5
    6  61,940,000-62,170,000      unmapped_end                   4
    6  61,980,000-62,370,000      split                         20
    6  95,500,000-95,740,000      unmapped_end                   5
    6  95,780,000-96,020,000      unmapped_start                 5
    6  107,110,000-107,500,000    length_out_of_tolerance       20
    6  157,370,000-157,600,000    unmapped_end                   4
    6  157,410,000-157,840,000    split                         24
    6  157,650,000-157,880,000    unmapped_start                 4
    6  167,760,000-168,000,000    unmapped_end                   5
    6  167,990,000-168,230,000    unmapped_start                 5
    6  170,090,000-170,320,000    unmapped_end                   4
    6  170,130,000-170,470,000    interior_deleted              15
    6  170,280,000-170,510,000    unmapped_start                 4
    6  170,870,000-171,110,000    unmapped_end                   5
    7  50,000-280,000             unmapped_end                   4
    7  90,000-430,000             interior_deleted              15
    7  240,000-470,000            unmapped_start                 4
    7  50,190,000-50,410,000      unmapped_end                   3
    7  50,220,000-50,570,000      interior_deleted              16
    7  50,380,000-50,600,000      unmapped_start                 3
    7  57,870,000-58,110,000      unmapped_end                   5
    7  61,000,000-61,240,000      unmapped_start                 5
    7  61,130,000-61,360,000      unmapped_end                   4
    7  61,170,000-61,460,000      interior_deleted              10
    7  61,270,000-61,510,000      unmapped_end                   5
    7  61,320,000-62,110,000      split                         60
    7  61,920,000-62,150,000      unmapped_start                 4
    7  72,300,000-72,690,000      length_out_of_tolerance       20
    7  74,530,000-75,110,000      split                         39
    7  98,190,000-98,580,000      length_out_of_tolerance       20
    7  100,370,000-100,600,000    unmapped_end                   4
    7  100,410,000-100,800,000    split                         20
    7  129,970,000-130,210,000    unmapped_end                   5
    7  130,200,000-130,440,000    split                          5
    7  139,190,000-139,400,000    unmapped_end                   2
    7  139,210,000-139,570,000    interior_deleted              17
    7  139,380,000-139,590,000    unmapped_start                 2
    7  141,860,000-142,090,000    unmapped_end                   4
    7  141,900,000-142,240,000    split                         15
    7  142,050,000-142,270,000    unmapped_start                 3
    7  142,080,000-142,290,000    unmapped_both                  2
    7  142,100,000-142,320,000    unmapped_end                   3
    7  142,130,000-142,670,000    split                         35
    7  143,160,000-143,390,000    unmapped_end                   4
    7  143,200,000-143,540,000    interior_deleted              15
    7  143,350,000-143,580,000    unmapped_start                 4
    7  154,090,000-154,330,000    unmapped_end                   5
    7  154,320,000-154,560,000    unmapped_start                 5
    8  2,100,000-2,510,000        split                         22
    8  2,320,000-2,520,000        length_out_of_tolerance        1
    8  7,290,000-7,520,000        unmapped_end                   4
    8  7,330,000-7,670,000        interior_deleted              15
    8  7,480,000-7,710,000        unmapped_start                 4
    8  11,910,000-12,140,000      unmapped_end                   4
    8  11,950,000-12,290,000      interior_deleted              15
    8  12,100,000-12,330,000      unmapped_start                 4
    8  43,650,000-43,890,000      unmapped_end                   5
    8  46,780,000-47,020,000      unmapped_start                 5
    8  48,030,000-48,410,000      split                         19
    8  86,390,000-86,630,000      unmapped_end                   5
    8  86,670,000-86,910,000      unmapped_start                 5
    8  142,560,000-143,000,000    split                         25
    8  144,900,000-145,130,000    length_out_of_tolerance        4
    8  145,140,000-145,680,000    split                         17
    8  146,120,000-146,360,000    unmapped_end                   5
    9  39,260,000-40,660,000      split                        117
    9  40,750,000-41,560,000      split                         62
    9  41,370,000-41,590,000      unmapped_start                 3
    9  41,400,000-41,780,000      split                         19
    9  42,430,000-42,660,000      unmapped_end                   4
    9  42,470,000-42,860,000      split                         20
    9  42,670,000-43,030,000      length_out_of_tolerance       17
    9  42,860,000-44,140,000      split                         84
    9  43,950,000-44,170,000      unmapped_start                 3
    9  43,980,000-44,870,000      split                         66
    9  44,680,000-44,900,000      unmapped_start                 3
    9  44,710,000-44,920,000      unmapped_both                  2
    9  44,730,000-44,950,000      unmapped_end                   3
    9  44,760,000-45,100,000      interior_deleted              15
    9  44,910,000-45,140,000      unmapped_start                 4
    9  45,070,000-45,310,000      unmapped_end                   5
    9  45,300,000-46,140,000      split                         65
    9  46,030,000-46,260,000      unmapped_end                   4
    9  46,070,000-46,520,000      split                         26
    9  46,510,000-46,750,000      unmapped_start                 5
    9  46,880,000-47,120,000      unmapped_end                   5
    9  47,110,000-47,370,000      split                          7
    9  65,410,000-66,110,000      split                         51
    9  65,920,000-66,150,000      unmapped_start                 4
    9  66,010,000-66,240,000      unmapped_end                   4
    9  66,050,000-66,600,000      split                         36
    9  66,410,000-66,610,000      unmapped_start                 1
    9  66,420,000-66,650,000      unmapped_both                  4
    9  66,460,000-66,660,000      unmapped_end                   1
    9  66,470,000-66,810,000      split                         15
    9  66,620,000-66,850,000      unmapped_start                 4
    9  66,680,000-66,910,000      unmapped_end                   4
    9  66,720,000-67,060,000      split                         15
    9  66,870,000-67,100,000      unmapped_start                 4
    9  66,920,000-67,160,000      unmapped_end                   5
    9  67,150,000-67,360,000      unmapped_start                 2
    9  67,170,000-67,400,000      unmapped_both                  4
    9  67,210,000-67,420,000      unmapped_end                   2
    9  67,460,000-68,040,000      split                         39
    9  68,080,000-68,320,000      unmapped_start                 5
    9  68,330,000-68,570,000      unmapped_end                   5
    9  68,610,000-68,890,000      split                          9
    9  68,930,000-69,170,000      unmapped_start                 5
    9  69,030,000-69,470,000      split                         25
    9  69,280,000-69,510,000      unmapped_start                 4
    9  69,830,000-70,060,000      unmapped_end                   4
    9  69,870,000-70,210,000      interior_deleted              15
    9  70,020,000-70,260,000      unmapped_both                  5
    9  70,070,000-70,270,000      unmapped_end                   1
    9  70,260,000-70,520,000      split                          7
    9  70,330,000-70,550,000      unmapped_end                   3
    9  70,360,000-70,790,000      split                         24
    9  70,780,000-71,020,000      unmapped_start                 5
    9  92,160,000-92,400,000      unmapped_end                   5
    9  92,620,000-92,860,000      unmapped_start                 5
    9  132,890,000-133,130,000    unmapped_end                   5
    9  133,170,000-133,410,000    unmapped_start                 5
    9  136,860,000-137,090,000    unmapped_end                   4
    9  136,900,000-137,240,000    interior_deleted              15
    9  137,050,000-137,280,000    unmapped_start                 4
    9  138,980,000-139,210,000    unmapped_end                   4
    9  139,020,000-139,360,000    interior_deleted              15
    9  139,170,000-139,400,000    unmapped_start                 4
    9  140,970,000-141,210,000    unmapped_end                   5
   10  10,000-250,000             unmapped_start                 5
   10  17,700,000-18,310,000      split                         42
   10  38,630,000-38,860,000      unmapped_end                   4
   10  38,670,000-39,010,000      interior_deleted              15
   10  38,820,000-39,050,000      unmapped_start                 4
   10  38,970,000-39,210,000      unmapped_end                   5
   10  42,300,000-42,540,000      unmapped_start                 5
   10  42,360,000-42,590,000      unmapped_end                   4
   10  42,400,000-42,780,000      split                         19
   10  46,240,000-46,470,000      unmapped_end                   4
   10  46,280,000-46,660,000      split                         19
   10  46,840,000-47,480,000      split                         45
   10  47,470,000-47,710,000      unmapped_start                 5
   10  47,590,000-48,420,000      split                         46
   10  48,720,000-49,150,000      split                         24
   10  49,140,000-49,380,000      unmapped_start                 5
   10  50,950,000-51,590,000      split                         45
   10  51,400,000-51,630,000      unmapped_start                 4
   10  51,530,000-51,900,000      split                         18
   10  125,680,000-125,910,000    unmapped_end                   4
   10  125,720,000-126,060,000    interior_deleted              15
   10  125,870,000-126,100,000    unmapped_start                 4
   10  128,430,000-128,670,000    unmapped_end                   5
   10  128,710,000-128,950,000    unmapped_start                 5
   10  133,200,000-133,430,000    unmapped_end                   4
   10  133,240,000-133,580,000    interior_deleted              15
   10  133,390,000-133,620,000    unmapped_start                 4
   10  133,490,000-133,720,000    unmapped_end                   4
   10  133,530,000-133,870,000    interior_deleted              15
   10  133,680,000-133,910,000    unmapped_start                 4
   10  135,280,000-135,530,000    length_out_of_tolerance        6
   11  10,000-240,000             unmapped_start                 4
   11  980,000-1,350,000          split                         18
   11  1,160,000-1,360,000        interior_deleted               1
   11  1,170,000-1,400,000        unmapped_start                 4
   11  49,840,000-50,230,000      length_out_of_tolerance       20
   11  50,600,000-50,840,000      split                          5
   11  51,040,000-51,280,000      split                          5
   11  51,410,000-51,650,000      unmapped_end                   5
   11  54,640,000-54,880,000      unmapped_start                 5
   11  68,900,000-69,130,000      unmapped_end                   4
   11  68,940,000-69,280,000      interior_deleted              15
   11  69,090,000-69,320,000      unmapped_start                 4
   11  69,540,000-69,770,000      unmapped_end                   4
   11  69,580,000-69,920,000      interior_deleted              15
   11  69,730,000-69,960,000      unmapped_start                 4
   11  70,610,000-71,000,000      length_out_of_tolerance       20
   11  87,500,000-87,880,000      split                         19
   11  87,690,000-87,920,000      unmapped_start                 4
   11  96,100,000-96,340,000      unmapped_end                   5
   11  96,380,000-96,620,000      unmapped_start                 5
   11  134,760,000-135,000,000    unmapped_end                   5
   12  10,000-340,000             split                         14
   12  7,000,000-7,230,000        unmapped_end                   4
   12  7,040,000-7,380,000        interior_deleted              15
   12  7,190,000-7,420,000        unmapped_start                 4
   12  34,670,000-34,910,000      unmapped_end                   5
   12  37,800,000-38,040,000      unmapped_start                 5
   12  109,190,000-109,420,000    unmapped_end                   4
   12  109,230,000-109,570,000    interior_deleted              15
   12  109,380,000-109,610,000    unmapped_start                 4
   12  122,350,000-122,580,000    unmapped_end                   4
   12  122,390,000-122,730,000    interior_deleted              15
   12  122,540,000-122,770,000    unmapped_start                 4
   12  132,520,000-132,760,000    unmapped_end                   5
   12  132,750,000-133,000,000    split                          6
   13  18,970,000-19,210,000      unmapped_start                 5
   13  86,580,000-86,820,000      unmapped_end                   5
   13  86,860,000-87,110,000      unmapped_start                 6
   13  112,170,000-112,410,000    unmapped_end                   5
   13  112,450,000-112,690,000    unmapped_start                 5
   13  114,140,000-114,380,000    unmapped_end                   5
   13  114,370,000-114,610,000    unmapped_start                 5
   13  114,450,000-114,690,000    unmapped_end                   5
   13  114,680,000-114,920,000    unmapped_start                 5
   13  114,920,000-115,160,000    unmapped_end                   5
   14  18,950,000-19,180,000      unmapped_start                 4
   14  19,290,000-19,680,000      length_out_of_tolerance       20
   14  19,570,000-20,010,000      split                         25
   14  19,900,000-20,290,000      length_out_of_tolerance       20
   15  19,950,000-20,180,000      unmapped_start                 4
   15  20,710,000-21,010,000      split                         11
   15  20,820,000-21,090,000      interior_deleted               8
   15  20,900,000-21,120,000      unmapped_start                 3
   15  21,210,000-21,450,000      unmapped_end                   5
   15  21,830,000-22,040,000      unmapped_start                 2
   15  21,850,000-22,070,000      split                          3
   15  23,220,000-23,750,000      split                         34
   15  28,970,000-29,200,000      unmapped_end                   4
   15  29,010,000-29,350,000      interior_deleted              15
   15  29,160,000-29,390,000      unmapped_start                 4
   15  82,450,000-83,220,000      split                         58
   15  84,740,000-84,990,000      length_out_of_tolerance        6
   15  84,800,000-85,030,000      unmapped_end                   4
   15  84,840,000-85,220,000      split                         19
   16  10,000-250,000             unmapped_start                 5
   16  8,450,000-8,680,000        unmapped_end                   4
   16  8,490,000-8,830,000        interior_deleted              15
   16  8,640,000-8,870,000        unmapped_start                 4
   16  18,340,000-18,730,000      length_out_of_tolerance       20
   16  33,150,000-33,490,000      length_out_of_tolerance       15
   16  33,840,000-34,080,000      unmapped_end                   5
   16  34,120,000-34,360,000      unmapped_start                 5
   16  35,100,000-35,340,000      unmapped_end                   5
   16  46,330,000-46,600,000      split                          8
   16  88,200,000-88,430,000      unmapped_end                   4
   16  88,240,000-88,580,000      interior_deleted              15
   16  88,390,000-88,620,000      unmapped_start                 4
   16  90,110,000-90,350,000      unmapped_end                   5
   17  110,000-350,000            unmapped_end                   5
   17  340,000-580,000            unmapped_start                 5
   17  21,350,000-21,620,000      split                          8
   17  21,610,000-21,850,000      unmapped_start                 5
   17  21,710,000-22,100,000      length_out_of_tolerance       20
   17  22,080,000-22,320,000      unmapped_end                   5
   17  25,210,000-25,450,000      unmapped_start                 5
   17  34,400,000-34,910,000      split                         32
   17  36,160,000-36,550,000      length_out_of_tolerance       20
   17  41,190,000-41,580,000      length_out_of_tolerance       20
   17  62,230,000-62,460,000      unmapped_end                   4
   17  62,270,000-62,610,000      interior_deleted              15
   17  62,420,000-62,650,000      unmapped_start                 4
   17  77,360,000-77,590,000      unmapped_end                   4
   17  77,400,000-77,820,000      split                         23
   17  79,520,000-79,750,000      unmapped_end                   4
   17  79,560,000-79,900,000      interior_deleted              15
   17  79,710,000-79,940,000      unmapped_start                 4
   18  15,230,000-15,470,000      unmapped_end                   5
   18  18,460,000-18,700,000      unmapped_start                 5
   18  44,350,000-44,740,000      length_out_of_tolerance       20
   18  51,870,000-52,110,000      unmapped_end                   5
   18  52,150,000-52,390,000      unmapped_start                 5
   18  72,100,000-72,330,000      unmapped_end                   4
   18  72,140,000-72,480,000      interior_deleted              15
   18  72,290,000-72,520,000      unmapped_start                 4
   18  75,540,000-75,770,000      unmapped_end                   4
   18  75,580,000-75,920,000      interior_deleted              15
   18  75,730,000-75,960,000      unmapped_start                 4
   18  77,830,000-78,070,000      split                          5
   19  10,000-240,000             unmapped_start                 4
   19  7,120,000-7,540,000        split                         23
   19  7,350,000-7,580,000        unmapped_start                 4
   19  8,500,000-8,730,000        unmapped_end                   4
   19  8,540,000-8,770,000        interior_deleted               4
   19  8,580,000-8,930,000        split                         16
   19  20,340,000-20,570,000      unmapped_end                   4
   19  20,380,000-20,780,000      split                         21
   19  24,450,000-24,690,000      unmapped_end                   5
   19  27,680,000-27,920,000      unmapped_start                 5
   19  40,210,000-40,600,000      split                         20
   20  10,000-240,000             unmapped_start                 4
   20  26,130,000-26,370,000      unmapped_end                   5
   20  29,360,000-29,600,000      unmapped_start                 5
   20  29,470,000-29,710,000      unmapped_end                   5
   20  29,750,000-29,990,000      unmapped_start                 5
   20  34,710,000-34,940,000      unmapped_end                   4
   20  34,750,000-35,090,000      interior_deleted              15
   20  34,900,000-35,130,000      unmapped_start                 4
   20  53,930,000-54,320,000      length_out_of_tolerance       20
   20  60,910,000-61,140,000      unmapped_end                   4
   20  60,950,000-61,210,000      interior_deleted               7
   20  61,020,000-61,260,000      unmapped_end                   5
   20  61,070,000-61,290,000      interior_deleted               3
   20  61,100,000-61,340,000      unmapped_start                 5
   20  61,150,000-61,410,000      interior_deleted               7
   20  61,220,000-61,450,000      unmapped_start                 4
   20  62,780,000-63,020,000      unmapped_end                   5
   21  9,360,000-9,590,000        unmapped_start                 4
   21  9,400,000-9,610,000        unmapped_both                  2
   21  9,420,000-9,640,000        unmapped_end                   3
   21  9,450,000-9,770,000        interior_deleted              13
   21  9,580,000-9,790,000        unmapped_end                   2
   21  9,600,000-9,820,000        unmapped_both                  3
   21  9,630,000-9,840,000        unmapped_start                 2
   21  9,650,000-9,970,000        interior_deleted              13
   21  9,780,000-10,010,000       unmapped_start                 4
   21  9,850,000-10,080,000       unmapped_end                   4
   21  9,890,000-10,210,000       interior_deleted              13
   21  10,020,000-10,230,000      unmapped_end                   2
   21  10,040,000-10,270,000      unmapped_both                  4
   21  10,310,000-10,550,000      unmapped_start                 5
   21  10,460,000-10,690,000      unmapped_end                   4
   21  10,500,000-10,840,000      split                         15
   21  10,650,000-10,880,000      unmapped_start                 4
   21  11,000,000-11,240,000      unmapped_end                   5
   21  14,280,000-14,520,000      unmapped_start                 5
   21  42,770,000-43,000,000      unmapped_end                   4
   21  42,810,000-43,150,000      interior_deleted              15
   21  42,960,000-43,190,000      unmapped_start                 4
   21  44,450,000-44,680,000      unmapped_end                   4
   21  44,490,000-44,830,000      interior_deleted              15
   21  44,640,000-44,870,000      unmapped_start                 4
   22  16,000,000-16,240,000      unmapped_start                 5
   22  16,510,000-16,750,000      unmapped_end                   5
   22  16,790,000-17,030,000      unmapped_start                 5
   22  18,500,000-18,890,000      length_out_of_tolerance       20
   22  20,150,000-20,530,000      split                         19
   22  20,340,000-20,560,000      unmapped_end                   3
   22  20,550,000-20,880,000      split                         14
   22  24,160,000-24,580,000      split                         23
   22  50,180,000-50,410,000      unmapped_end                   4
   22  50,220,000-50,560,000      interior_deleted              15
   22  50,370,000-50,600,000      unmapped_start                 4
   22  51,060,000-51,300,000      unmapped_end                   5

==============================================================================
DECISION: 1.826% of windows fail to lift cleanly vs the ~2% descope trigger -> windows product SHIPS (hpa1.12 proceeds)
==============================================================================
```
