scRNA-Seq Anlaysis - 01 - Processing
================

-   [Dependencies](#dependencies)
-   [Load the data](#load-the-data)
-   [Move data into scater](#move-data-into-scater)
    -   [colData](#coldata)
    -   [Setting up the SingleCellExperiment](#setting-up-the-singlecellexperiment)
    -   [Adding rowData](#adding-rowdata)
    -   [Adding exprs (log2+1) slot](#adding-exprs-log21-slot)
-   [Save the data](#save-the-data)

The following scRNA-Seq data was generated with the Fluidigm C1 HT 3' mRNA-Seq protocol. The samples are control and estrogen-treated primary cultures of ovarian epithelial cells. Fastq files were demultiplexed using Fluidigm's API, and the resulting fastq files for each individual capture site were processed with Kallisto.

Dependencies
============

``` r
library(scater)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colMeans,
    ##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax,
    ##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
    ##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
    ##     tapply, union, unique, unsplit, which, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: ggplot2

    ## Loading required package: SingleCellExperiment

    ## Loading required package: SummarizedExperiment

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:base':
    ## 
    ##     apply

    ## 
    ## Attaching package: 'scater'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(tximport)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:scater':
    ## 
    ##     arrange, filter, mutate, rename

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(EnsDb.Mmusculus.v79)
```

    ## Loading required package: ensembldb

    ## Loading required package: GenomicFeatures

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: AnnotationFilter

    ## 
    ## Attaching package: 'ensembldb'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     filter

    ## The following object is masked from 'package:scater':
    ## 
    ##     filter

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(biomaRt)
library(readr)
```

Load the data
=============

Getting the file names ready

``` r
row.id <- paste0("ROW", 1:40)
column.id <- rep(1:20, 40) %>% sort()
zero.padding <- c(rep("0", 40*9), rep("", 40*11))
column.id <- paste0(zero.padding, column.id)
column.id<- paste0("COL", column.id, "_")
samples <- paste0(column.id, row.id)
directories <- samples
head(samples)
```

    ## [1] "COL01_ROW1" "COL01_ROW2" "COL01_ROW3" "COL01_ROW4" "COL01_ROW5"
    ## [6] "COL01_ROW6"

Loading it in

``` r
edb <- EnsDb.Mmusculus.v79
Tx <- transcripts(edb, return.type="DataFrame")
tx2gene <- Tx[,c(1,7)]

files <- file.path(samples, "abundance.tsv")
txi <- tximport(paste0("/Volumes/ExternalHD/2017-E2-scRNA-Seq/output/kallisto/",files),
                type="kallisto", tx2gene=tx2gene,
                countsFromAbundance="no")
```

    ## Note: importing `abundance.h5` is typically faster than `abundance.tsv`

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396 397 398 399 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421 422 423 424 425 426 427 428 429 430 431 432 433 434 435 436 437 438 439 440 441 442 443 444 445 446 447 448 449 450 451 452 453 454 455 456 457 458 459 460 461 462 463 464 465 466 467 468 469 470 471 472 473 474 475 476 477 478 479 480 481 482 483 484 485 486 487 488 489 490 491 492 493 494 495 496 497 498 499 500 501 502 503 504 505 506 507 508 509 510 511 512 513 514 515 516 517 518 519 520 521 522 523 524 525 526 527 528 529 530 531 532 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548 549 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569 570 571 572 573 574 575 576 577 578 579 580 581 582 583 584 585 586 587 588 589 590 591 592 593 594 595 596 597 598 599 600 601 602 603 604 605 606 607 608 609 610 611 612 613 614 615 616 617 618 619 620 621 622 623 624 625 626 627 628 629 630 631 632 633 634 635 636 637 638 639 640 641 642 643 644 645 646 647 648 649 650 651 652 653 654 655 656 657 658 659 660 661 662 663 664 665 666 667 668 669 670 671 672 673 674 675 676 677 678 679 680 681 682 683 684 685 686 687 688 689 690 691 692 693 694 695 696 697 698 699 700 701 702 703 704 705 706 707 708 709 710 711 712 713 714 715 716 717 718 719 720 721 722 723 724 725 726 727 728 729 730 731 732 733 734 735 736 737 738 739 740 741 742 743 744 745 746 747 748 749 750 751 752 753 754 755 756 757 758 759 760 761 762 763 764 765 766 767 768 769 770 771 772 773 774 775 776 777 778 779 780 781 782 783 784 785 786 787 788 789 790 791 792 793 794 795 796 797 798 799 800 
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
colnames(txi$counts) <- samples
colnames(txi$abundance) <- samples
rm(tx2gene, files, row.id, column.id, zero.padding, directories)
```

Move data into scater
=====================

Scater takes the expression matrix, colData (cell attributes), and rowData (gene attributes)

colData
-------

``` r
cell.info <- as.data.frame(samples)
rownames(cell.info) <- samples
colnames(cell.info) <- "Cell"
cell.info$Condition <- c(rep("Control", 400), rep("Estrogen", 400))
cell.info$Condition <- factor(cell.info$Condition, levels=unique(cell.info$Condition))
cell.info$IFC.Column <- rep(1:20, 40) %>% sort(decreasing=F)
cell.info$IFC.Column <- factor(cell.info$IFC.Column, levels=unique(cell.info$IFC.Column))
cell.info$IFC.Row <- rep(1:40, 20)
cell.info$IFC.Row <- factor(cell.info$IFC.Row, levels=unique(cell.info$IFC.Row))
rownames(cell.info) <- cell.info$Cell
#Cell number and visual quality assessment of each cell was manually entered into an Excel file
#during microscopic inspection of the IFC after capture
ifc.annotation <- read.table("../data/IFC.Annotation.txt",
                             sep="\t", header=T)
cell.info$CellNumber <- ifc.annotation$CellNumber
```

Setting up the SingleCellExperiment
-----------------------------------

``` r
sce <- SingleCellExperiment(assays = list(counts=txi$counts),
                            colData=cell.info)
sce <- calculateQCMetrics(sce) 
```

Adding rowData
--------------

``` r
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
ensembl_symbols <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
                         filters="ensembl_gene_id",
                         values=rownames(counts(sce)),
                         mart=ensembl)
```

    ## Warning in strptime(x, fmt, tz = "GMT"): unknown timezone 'zone/tz/2017c.
    ## 1.0/zoneinfo/America/Toronto'

Annotating the genes a little

``` r
is.mito <- ensembl_symbols[grepl("^mt-",
                                 ensembl_symbols$mgi_symbol),]$ensembl_gene_id
is.protein <- dplyr::filter(as.data.frame(Tx), tx_biotype=="protein_coding")$gene_id

rowData(sce)$ensembl_gene_id <- rownames(counts(sce))
rowData(sce)$is.mito <- rownames(counts(sce)) %in% is.mito
rowData(sce)$is.protein <- rownames(counts(sce)) %in% is.protein

sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
```

Adding exprs (log2+1) slot
--------------------------

``` r
exprs(sce) <- log2(counts(sce) + 1)
```

Save the data
=============

``` r
saveRDS(sce, file="../output/sce.processed.rds")
```
