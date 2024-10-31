# Signal Analysis and Filtering Project

This MATLAB project processes and analyzes a pulse signal from two LEDs—Red and Infrared—recorded in a file, `pulse.txt`. The analysis includes signal loading, preprocessing, low-pass and high-pass filtering, pulse rate detection, and oxygen saturation (SaO₂) estimation.

## Project Structure

### 1. Parameters
The code defines key parameters:
- **Sampling frequency (`Fs`)**: 100 Hz.
- **Max pulse rate (`max_pulse_rate_hz`)**: Filters out high frequencies above 3 Hz (equivalent to 180 bpm).
- **Min breath rate (`min_breath_rate_hz`)**: Removes frequencies below 0.5 Hz to retain essential information in the Red signal.

### 2. Data Loading
The data is loaded from `pulse.txt` using:
- **`Led_R`**: Red LED signal.
- **`Led_IR`**: Infrared LED signal.

### 3. Preprocessing
- **Data Trimming**: Skips the first 10 seconds and then selects a 30-second observation window.

### 4. Signal Filtering
- **Low-pass Filtering**: Uses the FFT to remove frequencies above 3 Hz, isolating the pulse rate band.
- **High-pass Filtering**: Filters out frequencies below 0.5 Hz in the Red signal, preserving breath rate frequencies.

### 5. Peak Detection and Pulse Rate Estimation
- Peaks are detected in the high-pass filtered Red signal, and the average interval between peaks is used to calculate the pulse rate in BPM.

### 6. Oxygen Saturation (SaO₂) Estimation
- SaO₂ is estimated based on the AC and DC components of both Red and Infrared signals using interpolation.

### 7. Pulse Rate Estimation Using FFT
- An alternative pulse rate estimation method uses the FFT of the filtered Red signal in the 60–90 BPM range.

## Plots and Output

- **Original vs Filtered Signals**: Compares original and low-pass filtered signals for both Red and Infrared LEDs.
- **High-pass Filtered Signal with Peaks**: Shows the high-pass filtered Red signal with detected peaks.
- **Interpolated Max/Min Curves**: Displays interpolated curves for estimating SaO₂.
- **Pulse Rate via FFT**: Highlights the peak in the FFT corresponding to pulse rate.

## Example Output

The script outputs the following values:
- **Estimated Pulse Rate (BPM)**
- **Estimated SaO₂ (%)**

## Dependencies
- MATLAB or MATLAB-compatible environment.
- `pulse.txt` data file in the same directory as the code.

## Running the Code
Clone this repository, place `pulse.txt` in the root directory, and execute the code in MATLAB:
```matlab
run('your_script_name.m')
