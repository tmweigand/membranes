# Instructions to perform Coordination number analysis

## Example commands

```bash
python coordination_number/coordination_number.py --system vickers_membrane --input data/water_data/ --output data_out/vickers_membrane/

python coordination_number/coordination_number.py --system xu_opls_9_ss --input <path-to-xu-data>/ --output data_out/xu_opls_9_ss/
```

## Arguments

1. `--system`: the system to analyze. Choose one of:

   | System | Description |
   |---|---|
   | `water_box` | Water box |
   | `vickers_membrane` | Vickers membrane |
   | `xu_cfvv_9_ss` | Xu et al. 2026, CFVV forcefield, 9 nm domain, SS restriction |
   | `xu_cvff_20_ss` | Xu et al. 2026, CFVV forcefield, 20 nm domain, SS restriction |
   | `xu_gaff_9_ss` | Xu et al. 2026, GAFF forcefield, 9 nm domain, SS restriction |
   | `xu_opls_9_ss` | Xu et al. 2026, OPLS forcefield, 9 nm domain, SS restriction |

2. `--input`: directory containing the simulation data.
3. `--output`: directory for the results.

## Xu data

The Xu data comes from https://github.com/hxu97/TransportPressure/tree/main.
Clone the repository and point `--input` at the folder for the system you want.
The raw files are read as they are, with no manual editing.