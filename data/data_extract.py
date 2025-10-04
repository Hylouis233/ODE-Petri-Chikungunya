import cv2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import os
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

def load_image(image_path):
    image = cv2.imread(image_path)
    if image is None:
        raise ValueError(f"Cannot load image: {image_path}")
    return image

def detect_blue_red_regions(image):
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    
    lower_blue = np.array([100, 50, 50])
    upper_blue = np.array([130, 255, 255])
    blue_mask = cv2.inRange(hsv, lower_blue, upper_blue)
    
    lower_red1 = np.array([0, 50, 50])
    upper_red1 = np.array([10, 255, 255])
    lower_red2 = np.array([170, 50, 50])
    upper_red2 = np.array([180, 255, 255])
    
    red_mask1 = cv2.inRange(hsv, lower_red1, upper_red1)
    red_mask2 = cv2.inRange(hsv, lower_red2, upper_red2)
    red_mask = cv2.bitwise_or(red_mask1, red_mask2)
    
    return blue_mask, red_mask

def visualize_detection(image, blue_mask, red_mask):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    axes[0, 0].imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
    axes[0, 0].set_title('Original Image')
    axes[0, 0].axis('off')
    
    axes[0, 1].imshow(blue_mask, cmap='gray')
    axes[0, 1].set_title('Blue Region Detection')
    axes[0, 1].axis('off')
    
    axes[0, 2].imshow(red_mask, cmap='gray')
    axes[0, 2].set_title('Red Region Detection')
    axes[0, 2].axis('off')
    
    blue_overlay = image.copy()
    blue_overlay[blue_mask > 0] = [255, 0, 0]
    axes[1, 0].imshow(cv2.cvtColor(blue_overlay, cv2.COLOR_BGR2RGB))
    axes[1, 0].set_title('Blue Region Overlay')
    axes[1, 0].axis('off')
    
    red_overlay = image.copy()
    red_overlay[red_mask > 0] = [0, 255, 0]
    axes[1, 1].imshow(cv2.cvtColor(red_overlay, cv2.COLOR_BGR2RGB))
    axes[1, 1].set_title('Red Region Overlay')
    axes[1, 1].axis('off')
    
    combined_overlay = image.copy()
    combined_overlay[blue_mask > 0] = [255, 0, 0]
    combined_overlay[red_mask > 0] = [0, 255, 0]
    axes[1, 2].imshow(cv2.cvtColor(combined_overlay, cv2.COLOR_BGR2RGB))
    axes[1, 2].set_title('Combined Detection Result')
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    plt.show()

def export_bars_to_csv(blue_bars, red_bars, image_name):
    if blue_bars:
        blue_data = []
        for i, bar in enumerate(blue_bars):
            blue_data.append({
                'blue_width': bar['width'],
                'blue_height': bar['height'],
                'blue_x_coordinate': bar['x'],
                'blue_area': bar['area']
            })
        
        blue_df = pd.DataFrame(blue_data)
        blue_filename = f'{image_name}_blue_bars_dimensions.csv'
        blue_df.to_csv(blue_filename, index=False, encoding='utf-8-sig')
    
    if red_bars:
        red_data = []
        for i, bar in enumerate(red_bars):
            red_data.append({
                'red_width': bar['width'],
                'red_height': bar['height'],
                'red_x_coordinate': bar['x'],
                'red_area': bar['area']
            })
        
        red_df = pd.DataFrame(red_data)
        red_filename = f'{image_name}_red_bars_dimensions.csv'
        red_df.to_csv(red_filename, index=False, encoding='utf-8-sig')
    
    if blue_bars and red_bars:
        comparison_data = []
        max_bars = max(len(blue_bars), len(red_bars))
        
        for i in range(max_bars):
            row = {}
            
            if i < len(blue_bars):
                row.update({
                    'blue_width': blue_bars[i]['width'],
                    'blue_height': blue_bars[i]['height'],
                    'blue_x_coordinate': blue_bars[i]['x'],
                    'blue_area': blue_bars[i]['area']
                })
            else:
                row.update({
                    'blue_width': None,
                    'blue_height': None,
                    'blue_x_coordinate': None,
                    'blue_area': None
                })
            
            if i < len(red_bars):
                row.update({
                    'red_width': red_bars[i]['width'],
                    'red_height': red_bars[i]['height'],
                    'red_x_coordinate': red_bars[i]['x'],
                    'red_area': red_bars[i]['area']
                })
            else:
                row.update({
                    'red_width': None,
                    'red_height': None,
                    'red_x_coordinate': None,
                    'red_area': None
                })
            
            comparison_data.append(row)
        
        comparison_df = pd.DataFrame(comparison_data)
        comparison_filename = f'{image_name}_bars_comparison.csv'
        comparison_df.to_csv(comparison_filename, index=False, encoding='utf-8-sig')

def save_masks_separately(blue_mask, red_mask, image_name):
    blue_mask_filename = f'{image_name}_blue_mask.png'
    cv2.imwrite(blue_mask_filename, blue_mask)
    
    red_mask_filename = f'{image_name}_red_mask.png'
    cv2.imwrite(red_mask_filename, red_mask)

def find_bar_contours(mask, color_name):
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    bars_info = []
    for i, contour in enumerate(contours):
        x, y, w, h = cv2.boundingRect(contour)
        
        if w > 5 and h > 5:
            bars_info.append({
                'id': i,
                'x': x,
                'y': y,
                'width': w,
                'height': h,
                'center_x': x + w // 2,
                'center_y': y + h // 2,
                'area': cv2.contourArea(contour)
            })
    
    bars_info.sort(key=lambda bar: bar['x'])
    
    return bars_info, contours

def create_annotated_image(image, blue_mask, red_mask, blue_bars, red_bars, blue_contours, red_contours):
    annotated = image.copy()
    
    cv2.drawContours(annotated, blue_contours, -1, (255, 255, 0), 2)
    for bar in blue_bars:
        cv2.rectangle(annotated, (bar['x'], bar['y']), 
                     (bar['x'] + bar['width'], bar['y'] + bar['height']), 
                     (255, 255, 0), 1)
        
        text = f"W:{bar['width']} H:{bar['height']}"
        text_x = bar['x']
        text_y = bar['y'] + bar['height'] + 15
        cv2.putText(annotated, text, (text_x, text_y), 
                   cv2.FONT_HERSHEY_SIMPLEX, 0.4, (255, 255, 0), 1)
    
    cv2.drawContours(annotated, red_contours, -1, (0, 255, 255), 2)
    for bar in red_bars:
        cv2.rectangle(annotated, (bar['x'], bar['y']), 
                     (bar['x'] + bar['width'], bar['y'] + bar['height']), 
                     (0, 255, 255), 1)
        
        text = f"W:{bar['width']} H:{bar['height']}"
        text_x = bar['x']
        text_y = bar['y'] + bar['height'] + 30
        cv2.putText(annotated, text, (text_x, text_y), 
                   cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 255, 255), 1)
    
    return annotated

def visualize_with_annotations(image, blue_mask, red_mask, blue_bars, red_bars, blue_contours, red_contours, image_name):
    annotated = create_annotated_image(image, blue_mask, red_mask, blue_bars, red_bars, blue_contours, red_contours)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    axes[0, 0].imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
    axes[0, 0].set_title('Original Image')
    axes[0, 0].axis('off')
    
    axes[0, 1].imshow(blue_mask, cmap='Blues')
    axes[0, 1].set_title('Blue Region Detection')
    axes[0, 1].axis('off')
    
    axes[1, 0].imshow(red_mask, cmap='Reds')
    axes[1, 0].set_title('Red Region Detection')
    axes[1, 0].axis('off')
    
    axes[1, 1].imshow(cv2.cvtColor(annotated, cv2.COLOR_BGR2RGB))
    axes[1, 1].set_title('Size Annotation Result')
    axes[1, 1].axis('off')
    
    plt.tight_layout()
    result_filename = f'{image_name}_detection_results_with_annotations.png'
    plt.savefig(result_filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    annotated_filename = f'{image_name}_annotated_chart.png'
    cv2.imwrite(annotated_filename, annotated)

def analyze_chart_structure(image, blue_mask, red_mask):
    height, width = image.shape[:2]
    
    blue_pixels = np.column_stack(np.where(blue_mask > 0))
    red_pixels = np.column_stack(np.where(red_mask > 0))
    
    if len(blue_pixels) > 0:
        blue_bounds = {
            'min_y': blue_pixels[:, 0].min(),
            'max_y': blue_pixels[:, 0].max(),
            'min_x': blue_pixels[:, 1].min(),
            'max_x': blue_pixels[:, 1].max()
        }
    
    if len(red_pixels) > 0:
        red_bounds = {
            'min_y': red_pixels[:, 0].min(),
            'max_y': red_pixels[:, 0].max(),
            'min_x': red_pixels[:, 1].min(),
            'max_x': red_pixels[:, 1].max()
        }

def process_single_image(image_path, image_name):
    image_for_analysis = load_image(image_path)
    blue_mask, red_mask = detect_blue_red_regions(image_for_analysis)
    save_masks_separately(blue_mask, red_mask, image_name)
    blue_bars, blue_contours = find_bar_contours(blue_mask, "blue")
    red_bars, red_contours = find_bar_contours(red_mask, "red")
    export_bars_to_csv(blue_bars, red_bars, image_name)
    visualize_detection(image_for_analysis, blue_mask, red_mask)
    visualize_with_annotations(image_for_analysis, blue_mask, red_mask, blue_bars, red_bars, blue_contours, red_contours, image_name)
    analyze_chart_structure(image_for_analysis, blue_mask, red_mask)
    return True

def main():
    image_files = [
        ("foshan_left.png", "foshan_left"),
        ("foshan_right.png", "foshan_right"),
        ("foshan.webp", "foshan_original")
    ]
    
    for image_path, image_name in image_files:
        if os.path.exists(image_path):
            process_single_image(image_path, image_name)

if __name__ == "__main__":
    main()