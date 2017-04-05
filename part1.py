import cv2
import numpy as np

print(cv2.__version__)

width = 100
height = 100

image = np.array([[0, 1, 0]])

cv2.imwrite('/tmp/image.ppm', image)
