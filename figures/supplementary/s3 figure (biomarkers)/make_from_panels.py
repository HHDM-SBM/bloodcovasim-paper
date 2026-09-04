from PIL import Image, ImageDraw, ImageFont
import numpy as np
import os

test_counts = [10, 50]
image_paths = [
    f'../../biomarkers_detection_rates_{count}.png' for count in test_counts
]

dpi = 300
width = 18
height = 16
width_pixels = int(width * dpi)
height_pixels = int(height * dpi)
panel_width = width_pixels
panel_height = height_pixels // 2
scale = 1.0

result = Image.new('RGBA', (width_pixels, height_pixels), (255, 255, 255, 255))
draw = ImageDraw.Draw(result)

img = Image.open(image_paths[0])
# w, h = height_pixels * img.size[0] // img.size[1], height_pixels 
w, h = panel_width, panel_width * img.size[1] // img.size[0]
w, h = int(w * scale), int(h * scale)
img = img.resize((w, h), Image.Resampling.LANCZOS)
result.paste(img, ((panel_width - w) // 2, (panel_height - h) // 2))


img = Image.open(image_paths[1])
# w, h = height_pixels * img.size[0] // img.size[1], height_pixels
w, h = panel_width, panel_width * img.size[1] // img.size[0]
w, h = int(w * scale), int(h * scale)
img = img.resize((w, h), Image.Resampling.LANCZOS)
result.paste(img, ((panel_width - w) // 2, panel_height + (panel_height - h) // 2))


padding = 150

try:
    font = ImageFont.truetype("C:\\Users\\Lenovo\\AppData\\Local\\Microsoft\\Windows\\Fonts\\Inter-VariableFont_opsz,wght.ttf", 125)
except:
    font = ImageFont.load_default()

panel_labels = ['A', 'B', 'C', 'D']

positions = [
    (padding, padding // 2),
    (panel_width // 2 + padding, padding // 2),
    (padding, panel_height + padding // 2),
    (panel_width // 2 + padding, panel_height + padding // 2),
]

for label, pos in zip(panel_labels, positions):
    draw.text(pos, label, fill='black', font=font, anchor='mm')


panel_titles = [
    f'{test} ({count} daily tests)' for count in test_counts for test in ['AIC', 'Anderson-Darling'] 
]

positions = [
    (panel_width // 4, padding // 2),
    (3 * panel_width // 4, padding // 2),
    (panel_width // 4, panel_height + padding // 2),
    (3 * panel_width // 4, panel_height + padding // 2),
]

for title, pos in zip(panel_titles, positions):
    draw.text(pos, title, fill='black', font=font, anchor='mm')

result.save('S3 Figure.png', dpi=(dpi, dpi))