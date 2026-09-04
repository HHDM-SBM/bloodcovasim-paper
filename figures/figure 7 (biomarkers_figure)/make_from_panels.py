from PIL import Image, ImageDraw, ImageFont
import numpy as np
import os

test_counts = [250]
image_paths = [
    f'../biomarkers_detection_rates_{count}.png' for count in test_counts
]

dpi = 300
width = 18
height = 8.5
width_pixels = int(width * dpi)
height_pixels = int(height * dpi)
panel_width = width_pixels
panel_height = height_pixels

result = Image.new('RGBA', (width_pixels, height_pixels), (255, 255, 255, 255))
draw = ImageDraw.Draw(result)

img = Image.open(image_paths[0])
w, h = width_pixels, width_pixels * img.size[1] // img.size[0]
w, h = int(w), int(h)
img = img.resize((w, h), Image.Resampling.LANCZOS)
result.paste(img, ((panel_width - w) // 2, (panel_height - h)))


# img = Image.open(image_paths[1])
# w, h = height_pixels * img.size[0] // img.size[1], height_pixels 
# w, h = int(w * scale), int(h * scale)
# img = img.resize((w, h), Image.Resampling.LANCZOS)
# result.paste(img, (panel_width + (panel_width - w) // 2, (panel_height - h) // 2))


padding = 150

try:
    font = ImageFont.truetype("C:\\Users\\Lenovo\\AppData\\Local\\Microsoft\\Windows\\Fonts\\Inter-VariableFont_opsz,wght.ttf", 125)
except:
    font = ImageFont.load_default()

panel_labels = ['A', 'B']

positions = [
    (padding, (panel_height - h) // 2 + padding),
    (panel_width // 2 + padding, (panel_height - h) // 2 + padding),
]

for label, pos in zip(panel_labels, positions):
    draw.text(pos, label, fill='black', font=font, anchor='mm')



panel_labels = ['AIC', 'Anderson-Darling']

positions = [
    (panel_width // 4, (panel_height - h) // 2 + padding),
    (3 * panel_width // 4, (panel_height - h) // 2 + padding),
]

for label, pos in zip(panel_labels, positions):
    draw.text(pos, label, fill='black', font=font, anchor='mm')



panel_titles = [
    f'{count} daily blood tests' for count in test_counts
]

positions = [
    (panel_width // 2, padding),
]

for title, pos in zip(panel_titles, positions):
    draw.text(pos, title, fill='black', font=font, anchor='mm')

result.save('Figure 7.png', dpi=(dpi, dpi))