## Inputs

### VCF Format Support
`Verix` is designed to be flexible, supporting multiple ways of representing complex structural variants. 

Each VCF record in the input can represent a single breakpoint, a single breakend (BND), or a whole CSV. 
To recover CSVs from these diverse inputs, `Verix` relies on three distinct input formats specified by the user:

* **`multi_rec`**: The CSV is split across multiple VCF lines. 
These records are linked together by a shared parent ID `SVID` (or user-defined) in an `INFO` field (e.g., `SVID=123`). 
All records with the same ID are merged into one global SV object.

```angular2html
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
chr1    1883193 sv875_1 N       <CUT>   100     PASS    END=1886146;OP_TYPE=CUT;GRAMMAR=ABC->b;VSET=4;SVLEN=2954;SVID=sv875;SVTYPE=delINVdel;SYMBOL=A   GT      0|1
chr1    1886147 sv875_0 N       <INV>   100     PASS    END=1895367;OP_TYPE=INV;GRAMMAR=ABC->b;VSET=4;SVLEN=9221;SVID=sv875;SVTYPE=delINVdel;SYMBOL=B   GT      0|1
chr1    1895368 sv875_2 N       <CUT>   100     PASS    END=1904067;OP_TYPE=CUT;GRAMMAR=ABC->b;VSET=4;SVLEN=8700;SVID=sv875;SVTYPE=delINVdel;SYMBOL=C   GT      0|1
```

* **`single_rec`**: The entire CSV is described on one line. The `BKPS` (or a user-defined) `INFO` field contains a list of internal breakend positions.
the accepted formats for this field are `TAG1-Pos1-Pos2-Pos3...,TAG2-Pos4-Pos5...` or 
`TAG1-CHROM1:Pos1-CHROM2:Pos2-CHROM3:Pos3...,TAG2-CHROM1:Pos1...`. 
If the chromosomes are not specified in the `BKPS` field, all breakends are assumed to belong to the record's chromosome.
```angular2html
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  all
chr1    1214032 3       N       <CSV>   98      PASS    END=1221669;SVLEN=7637;SVTYPE=DEL+INV;SUPPORT=10;BKPS=DEL:7634-1214030-1221668,INV:2926-1214035-1216963 GT:DR:DV        ./.:22:10
chr1    1468702 5       N       <CSV>   97      PASS    END=1487574;SVLEN=18872;SVTYPE=DEL+INV;SUPPORT=10;BKPS=DEL:18831-1468737-1487570,INV:8855-1474629-1483484       GT:DR:DV        ./.:41:10
```

* **`default`**: Ignore the previously described `INFO` fields.

```angular2html
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
chr11   132760392   BND979    N   ]chr11:3356943]N   36.5    PASS    PRECISE;SVTYPE=BND;STRANDS=-+;MAPQ=36.5  GT:VAF:hVAF:DR:DV   0/1:1.00:1.00,0.00,0.00:0:3
chr11  3356943 BND980    N   N[chr11:132760392[  36.5    PASS    PRECISE;SVTYPE=BND;STRANDS=-+;MAPQ=36.5  GT:VAF:hVAF:DR:DV   0/1:1.00:1.00,0.00,0.00:0:3
```
  
### Optional Fields
Regardless of the format, the following fields are processed:
* **SV Type**: The `SVTYPE` (or user-defined) info field. Defaults to `unspecified` if missing.
* **CHROM2**: For records where the end breakend is on a different chromosome (inter-chromosomal). 
If not specified both breakends are assumed to belong to the record's chromosome.
* **TARGET** / **TARGET_CHROM**: Alternative way to specify a dispersion target position and chromosome.
`Verix` parser looks for those fields and create an additional breakend if it is present in the record.
* **ALT / REF**: If the record is a BND, `Verix` automatically leverages the standard breakend notation to pair breakends.
In this notation, mate BNDs are identified via the `ALT` field (e.g., `N[chr2:2222[`). 
* **END**: End position of a record.

---

## Outputs
`Verix` output VCFs follows the `single_rec` format, reporting breakends in the `BKPS` `INFO` field
as `SVID-CHROM1:Pos1-CHROM2:Pos2-CHROM3:Pos3...`.

```angular2html
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
chr1    9975680 sv2     N       <DEL>   .       .       END=9979621;CHROM2=chr1;SVTYPE=DEL;MATCH_LIST=105|complete|0|1|0|Contiguous|chr1:9979621_chr1:9979621|chr1:9975678_chr1:9975679;OPTIMAL_MATCH_MODE=complete;OPTIMAL_MATCH_TYPE=DEL;OPTIMAL_MATCH_ID=105;OPTIMAL_MATCH_SCORE=0-1-0-Contiguous;BND_MATCH=105:9975678,105:9979621  GT      0/1
chr1    109586402       sv892   N       <delINVdel>     .       .       END=109600644;CHROM2=chr1;SVTYPE=delINVdel;BKPS=sv892-chr1:109594562-chr1:109596655;MATCH_LIST=severus_83|complete|0|2|0|Contiguous|chr1:109594564_chr1:109594562|chr1:109596655_chr1:109596655|chr1:109600644_chr1:109600644|chr1:109586401_chr1:109586401;OPTIMAL_MATCH_MODE=complete;OPTIMAL_MATCH_TYPE=BND;OPTIMAL_MATCH_ID=severus_83;OPTIMAL_MATCH_SCORE=0-2-0-Contiguous;BND_MATCH=severus_83:109586401,severus_83:109594564,severus_83:109596655,severus_83:109600644  GT      1/1
```
