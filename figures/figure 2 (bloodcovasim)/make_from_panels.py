from PIL import Image, ImageDraw, ImageFont
import numpy as np
import os

test_counts = [10, 1250]
image_paths = [
    'panel A (crp_dynamics)/crp_dynamics.png',
    'panel B (flowchart)/flowchart.png',
]

dpi = 300
scale = 0.9

img = [Image.open(path) for path in image_paths]
result = Image.new('RGBA', (img[0].size[0] + img[1].size[0], img[0].size[1]), (255, 255, 255, 255))
draw = ImageDraw.Draw(result)

img = Image.open(image_paths[0])
original_size = img.size
img = img.resize((int(img.size[0] * scale), int(img.size[1] * scale)), Image.Resampling.LANCZOS)
result.paste(img, ((original_size[0] - img.size[0]) // 2, (original_size[1] - img.size[1]) // 2))

img1 = Image.open(image_paths[1])
original_size1 = img1.size
img1 = img1.resize((int(img1.size[0] * scale), int(img1.size[1] * scale)), Image.Resampling.LANCZOS)
result.paste(img1, (original_size[0] + (original_size1[0] - img1.size[0]) // 2, (original_size1[1] - img1.size[1]) // 2))

padding = 150

try:
    font = ImageFont.truetype("C:\\Users\\Lenovo\\AppData\\Local\\Microsoft\\Windows\\Fonts\\Inter-VariableFont_opsz,wght.ttf", 125)
except:
    font = ImageFont.load_default()

panel_labels = ['A', 'B']

positions = [
    (padding, padding),
    (original_size[0] + padding, padding),
]

for label, pos in zip(panel_labels, positions):
    draw.text(pos, label, fill='black', font=font, anchor='mm')

result.save('Figure 2.png', dpi=(dpi, dpi))