# Flyback Converter Simulation – State-Space Averaged Model (CCM)

This repository contains a complete MATLAB implementation of a flyback converter operating in **Continuous Conduction Mode (CCM)** using the **state-space averaging technique**. The model is based on the IEEE reference:

> A. S. Raj et al., "Modelling of Flyback Converter Using State-Space Averaging Technique," CONECCT 2015, DOI: [10.1109/CONECCT.2015.7383871](https://doi.org/10.1109/CONECCT.2015.7383871).

The project provides a unified framework for **large-signal simulation**, **small-signal modeling**, and **frequency-domain stability analysis** of the flyback converter.

---

## Features

-  Two-mode state-space modeling (Switch ON / Switch OFF)

-  Large-signal averaged model

-  Small-signal control-to-output and input-to-output transfer functions

-  PWM-based time-domain simulation using RK4 numerical integration

-  Steady-state performance metrics (voltage, ripple, current, power, efficiency)

-  Frequency-domain analysis:

  - Bode plots

  - Nyquist diagrams

  - Gain and phase margins

  - Automatic extraction of poles and zeros

-  High-resolution PWM discretization

---

## State Variables

| Variable | Description |

|----------|-------------|

| iₗₘ | Magnetizing inductance current |

| vC | Output capacitor voltage |

---

## Model Inputs

| Input | Description |

|-------|-------------|

| Vᵢₙ | Input voltage |

| V_d | Diode forward voltage |

| d | Duty cycle (control input for small-signal analysis) |

---

## Simulation Parameters (Default)

| Parameter | Symbol | Value |

|-----------|--------|-------|

| Input Voltage | Vᵢₙ | 24 V |

| Switching Frequency | f_sw | 100 kHz |

| Duty Cycle | d | 0.8 |

| Load Resistance | R | 10 Ω |

| Output Capacitance | C | 47 µF |

| Magnetizing Inductance | Lₘ | 100 µH |

| Turns Ratio | n | 2 |

| Switch Resistance | R_sw | 0.1 Ω |

| Capacitor ESR | R_c | 0.05 Ω |

| Diode Drop | V_d | 0.7 V |

---

## Output Quantities

- Output voltage and ripple

- Magnetizing current (average, peak, ripple)

- Output power

- Estimated efficiency

- Small-signal transfer functions:

  - **Gvd(s)** – Control-to-output

  - **Gvg(s)** – Input-to-output

---

## Visualization

The script automatically generates:

-  Time-domain waveforms

-  Zoomed switching waveforms

-  Phase-plane trajectory

-  Instantaneous power and efficiency

-  Bode and Nyquist plots with stability margins

---

## Usage

1. Open the MATLAB script.

2. Modify the electrical parameters if needed.

3. Run the script.

```matlab

% Example: Run the simulation

flyback_ccm_simulation

```

---

## Requirements

- MATLAB (R2018b or later recommended)

- Control System Toolbox (for transfer function analysis)

---

## Author

LAOUAR Ouassim

CentraleSupélec

---

## Reference

If you use this code in your research, please cite:

```bibtex

@inproceedings{raj2015modelling,

  title={Modelling of Flyback Converter Using State-Space Averaging Technique},

  author={Raj, A. S. and others},

  booktitle={2015 IEEE International Conference on Electronics, Computing and Communication Technologies (CONECCT)},

  year={2015},

  doi={10.1109/CONECCT.2015.7383871}

}


