from PIL import Image, ImageDraw, ImageFont
import numpy as np
import os

test_counts = [50, 250]
image_paths = [
    f'../../synthetic_detection_plot_{count}.png' for count in test_counts
]

dpi = 300
width = 24
height = 12
width_pixels = int(width * dpi)
height_pixels = int(height * dpi)
panel_width = width_pixels // 2
panel_height = height_pixels
scale = 0.85

result = Image.new('RGBA', (width_pixels, height_pixels), (255, 255, 255, 255))
draw = ImageDraw.Draw(result)

img = Image.open(image_paths[0])
w, h = height_pixels * img.size[0] // img.size[1], height_pixels 
w, h = int(w * scale), int(h * scale)
img = img.resize((w, h), Image.Resampling.LANCZOS)
result.paste(img, ((panel_width - w) // 2, (panel_height - h) // 2))


img = Image.open(image_paths[1])
w, h = height_pixels * img.size[0] // img.size[1], height_pixels 
w, h = int(w * scale), int(h * scale)
img = img.resize((w, h), Image.Resampling.LANCZOS)
result.paste(img, (panel_width + (panel_width - w) // 2, (panel_height - h) // 2))


padding = 150

try:
    font = ImageFont.truetype("C:\\Users\\Lenovo\\AppData\\Local\\Microsoft\\Windows\\Fonts\\Inter-VariableFont_opsz,wght.ttf", 125)
except:
    font = ImageFont.load_default()

panel_labels = ['A', 'B', 'C', 'D']

positions = [
    (padding, padding),
    (panel_width + padding, padding),
    (padding, panel_height // 2),
    (panel_width + padding, panel_height // 2),
]

for label, pos in zip(panel_labels, positions):
    draw.text(pos, label, fill='black', font=font, anchor='mm')


panel_titles = [
    f'{count} daily blood tests' for count in test_counts
]

positions = [
    (panel_width // 2, padding),
    (panel_width + panel_width // 2, padding),
]

for title, pos in zip(panel_titles, positions):
    draw.text(pos, title, fill='black', font=font, anchor='mm')

result.save('S1 Figure.png', dpi=(dpi, dpi))