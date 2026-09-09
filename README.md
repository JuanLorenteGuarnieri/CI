# Computational Imaging

Implementation of computational imaging techniques developed as part of the Computational Imaging course at the University of Zaragoza.

The repository contains four practical projects covering RAW image processing, computational deblurring, HDR reconstruction and transient 3D reconstruction.

## Projects

### 1. RAW Image Processing

End-to-end processing pipeline for RAW images, including sensor linearization, Bayer demosaicing, white balancing, denoising, color enhancement, tone reproduction and image compression.

**Techniques**
- RAW sensor linearization
- Bayer demosaicing
- White World, Gray World and manual white balancing
- Gaussian, median and mean filtering
- HSV-based color enhancement
- Exposure and gamma-based tone reproduction
- PNG/JPEG export
- GPU acceleration with MATLAB Parallel Computing Toolbox and CPU fallback

[View Lab 1 →](./lab1)

---

### 2. Computational Image Deblurring

Image restoration pipeline based on physically motivated defocus simulation and frequency-domain deconvolution.

The implementation supports multiple aperture models and compares different reconstruction methods, including Wiener deconvolution and Richardson–Lucy.

**Techniques**
- Defocus and Gaussian noise simulation
- Aperture-based PSF modelling
- Frequency-domain Wiener deconvolution
- Richardson–Lucy reconstruction
- `1/f` power-spectrum prior
- RGB channel-wise restoration
- Alternative aperture models

[View Lab 2 →](./lab2)

---

### 3. HDR Imaging and Tone Mapping

High Dynamic Range imaging pipeline based on multi-exposure image stacks.

The project estimates the camera response function using the Debevec method, reconstructs an HDR radiance map and applies different tone-mapping operators for visualization.

**Techniques**
- Multi-exposure image processing
- Camera response function estimation
- Debevec HDR reconstruction
- Weighted least-squares / linear-system formulation
- Radiance map reconstruction
- Reinhard tone mapping
- Durand tone mapping
- Bilateral filtering in log-intensity space

[View Lab 3 →](./lab3)

---

### 4. Transient 3D Reconstruction

Volumetric reconstruction of hidden scenes from transient measurements using back-projection.

The implementation supports both confocal and non-confocal measurements and incorporates attenuation, foreshortening and temporal filtering corrections.

**Techniques**
- Transient imaging
- Confocal and non-confocal measurements
- 3D voxel-grid reconstruction
- Time-of-flight based back-projection
- Attenuation correction
- Foreshortening correction
- Phasor-based temporal filtering
- Morlet wavelet filtering
- Reconstruction parameter sweeps
- Execution-time evaluation

[View Lab 4 →](./lab4)

---

## Technologies

- **MATLAB**
- Image Processing Toolbox
- Parallel Computing Toolbox
- FFT-based signal and image processing
- HDR imaging
- Image restoration
- Inverse problems
- Transient imaging
- Volumetric reconstruction

## Repository Structure

```text
CI/
├── lab1/   # RAW image processing
├── lab2/   # Computational deblurring
├── lab3/   # HDR imaging and tone mapping
└── lab4/   # Transient 3D reconstruction
```

## Authors

**Juan Lorente Guarnieri**  
**Hugo Mateo Trejo**
