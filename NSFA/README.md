# Non-stationary fluctuation analysis code

Last update: *13.10.2022*

## Code files

`NSFA_Selector.m` is a dedicated script that opens a file of interest and plots it before allowing the user to select the waves that are in a close enough region to eachother. Not for cherrypicking, more for excluding obviously weird ones. Note: baseline subtraction performed here

`peakerOGS.m` is a modded version of Andy Penn's peaker script that detects peaks in traces, updated to interact smoothly with ephysIO's updated features.

`avgtracesOGS.m` is a modded version of Andy Penn's peaker script that averages all of the above traces, updated to interact smoothly with ephysIO's updated features.

`nsfaOGS.m` is a modded version of Andy Penn's nsfa script that performs non stationary fluctuation analysis on the above output, updated to interact smoothly with ephysIO's updated features.

**Note:** More detailed descriptions of these scripts actions can be found in their respective descriptions.

## Analysis Protocol (OGS Na$_V$ NSFA)

1.  Use `NSFA_Selector.m` to clean up individual recordings

2.  Make note of the filepath and any important information (eg, Genotype, faults etc) in the master documents `NSFA_Master.xlsx` stored in relevant onedrive folder

3.  On each successive recording, perform the following steps

    1.  Run `peakerOGS.m` on the cleaned up trace that you're happy to include. Settings to note include;

        | UI Input Text           | Value |
        |-------------------------|-------|
        | Cut off frequency (kHz) | 1000  |

    2.  Run `avgtracesOGS.m` immediately after the above step. Setings to note include;

        | UI Input Text          | Value          |
        |------------------------|----------------|
        | Alignment              | start          |
        | Relative peak fraction | 0.95           |
        | Reference align        | steepest slope |

    3.  Run `nsfaOGS.m` immediately after the above step. Settings to note include;

        | UI Input Text         | Value                 |
        |-----------------------|-----------------------|
        | Correlation over time | Should be weak        |
        | peak from;            | relative to amplitude |
        | ensemble variance     | independent           |
        | bins                  | 14                    |
        | peak relationship     | 1                     |
