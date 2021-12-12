#from scipy.misc import imread
from PIL import Image
import numpy as np

sources = ['cao1.jpg', 'ksiwek1.jpg', 'zhang1.jpg']
images = []
for source in sources:
    image = np.asarray(Image.open('img/' + source))
    images.append(image)


mix = images[0] * 0.1 + images[1] * 0.2 + images[2] * 0.7
im = Image.fromarray(np.uint8(mix))
im.save('img/mix/mix1.jpg')

mix = images[0] * 0.3 + images[1] * 0.4 + images[2] * 0.3
im = Image.fromarray(np.uint8(mix))
im.save('img/mix/mix2.jpg')

mix = images[0] * 0.6 + images[1] * 0.2 + images[2] * 0.2
im = Image.fromarray(np.uint8(mix))
im.save('img/mix/mix3.jpg')