# Data · 数据 · Daten

The PIV data are not part of this repository. Put the course dataset here so that the folder looks like this, or point `PIV_DATA` to wherever it is:

PIV 数据不在仓库里。请把课程数据集放到这里，使目录结构如下；也可以用环境变量 `PIV_DATA` 指向它所在的位置：

Die PIV-Daten sind nicht Teil dieses Repositorys. Legen Sie den Kursdatensatz hier ab, sodass der Ordner wie unten aussieht, oder lassen Sie `PIV_DATA` auf seinen Speicherort zeigen:

```
data/
└── FM-Lab-WS2627-Group1/
    ├── Calibration/
    │   ├── Generated000000.T000.D000.P000.H000.LA.TIF
    │   └── Generated000000.T000.D000.P000.H000.LB.TIF
    ├── Data/
    │   ├── PIV_Data_Cyl.mat     with cylinder / 有圆柱 / mit Zylinder
    │   └── PIV_Data_Bulk.mat    empty channel / 空水槽 / leerer Kanal
    ├── Samples/
    │   └── cyl_0000xx.T000.D000.P000.H000.LA.TIF / .LB.TIF   (11 image pairs / 11 对图像 / 11 Bildpaare)
    └── Help.txt
```

```matlab
setenv('PIV_DATA', 'D:/path/to/FM-Lab-WS2627-Group1')
```

## Contents of the .mat files · .mat 文件内容 · Inhalt der .mat-Dateien

Each file holds one struct `pivdata` · 每个文件含一个结构体 `pivdata` · Jede Datei enthält eine Struktur `pivdata`:

| Field · 字段 · Feld | Size · 尺寸 · Größe | Meaning · 含义 · Bedeutung |
|---|---|---|
| `IAsize` | 1 | interrogation window, 32 px / 判读窗大小 / Auswertefenster |
| `Overlap` | 1 | overlap factor 2, i.e. 16 px vector spacing / 重叠系数，矢量间距 16 px / Überlappungsfaktor, Vektorabstand 16 px |
| `cam` | 1 × 2 | camera resolution [1040 1392] / 相机分辨率 / Kameraauflösung |
| `dt` | 1 | pulse separation, 0.003 s / 两次曝光间隔 / Pulsabstand |
| `dimX`, `dimY` | 1 | 86 × 64 vectors / 矢量个数 / Anzahl der Vektoren |
| `pixX`, `pixY` | 64 × 86 | vector positions in pixels; `pixY` points up / 矢量位置（像素），pixY 向上 / Vektorpositionen in Pixeln, `pixY` zeigt nach oben |
| `dispX`, `dispY` | 64 × 86 × N<sub>t</sub> | particle displacement in pixels, NaN where invalid; N<sub>t</sub> = 4069 (cylinder) and 690 (empty channel), although `Help.txt` says 6000 / 粒子位移（像素），无效处为 NaN；N<sub>t</sub> = 4069（有圆柱）和 690（空水槽），`Help.txt` 中写的 6000 不准 / Partikelverschiebung in Pixeln, NaN wo ungültig; N<sub>t</sub> = 4069 (Zylinder) und 690 (leerer Kanal), obwohl `Help.txt` 6000 angibt |

Processing done before delivery (from `Help.txt`): background subtraction, multi-pass interrogation 64 → 32 px with peak-height validation, median filter, no substitution of removed vectors. The magnification and the frame rate are not given and are determined by the scripts (`s1_calibration.m`, `s1b_ring.m`, `s6_flowrate_timing.m`).

数据在交付前已经做过（见 `Help.txt`）：减背景、64 → 32 px 多级判读并做峰高验证、中值滤波、剔除的矢量不替补。放大率和采样率没有给出，由脚本确定（`s1_calibration.m`、`s1b_ring.m`、`s6_flowrate_timing.m`）。

Vorverarbeitung vor der Übergabe (laut `Help.txt`): Hintergrundabzug, mehrstufige Auswertung 64 → 32 px mit Prüfung der Peakhöhe, Medianfilter, entfernte Vektoren werden nicht ersetzt. Abbildungsmaßstab und Bildrate sind nicht angegeben und werden von den Skripten bestimmt (`s1_calibration.m`, `s1b_ring.m`, `s6_flowrate_timing.m`).
