import cv2
import pathlib
import imutils
import easyocr
import numpy as np
from google.colab.patches import cv2_imshow

def deskew_croped_receipt(cropped_image):#Crop_Image

  """
  This function is taking cropped_image as input
  rotate it 36 times from (0 to 350)degree at interval of 10 degree.
  Finding maximum confidence(accuracy) of each word detected by easy OCR


  Parameters
  ---------
  Crop_Image: array
    A 3D array of Image,cropped_Image form other function output

  Returns
  -------

  array
    this array contain confidence_of_words_accuracy,count_of_words at different angles from 0 to 350 degree

  """

  confidence_and_angle_array = []

  for angle in range(0,360,10):

    # Storing image array at diffrent angle
    rotated = imutils.rotate_bound(cropped_image,angle) #Crop_Image

    reader = easyocr.Reader(['en'])
    # performing OCR on each receipt rotated at diffrent angles
    res = reader.readtext(rotated)
    
    #calculating no. of characters detected by our easy OCR
    characters_count = len(res)
    
    # above 130 OCR characters are useless in scanned copy of receipt,
    # hypertuning of useful_words_count parameter to fetch good words in a receipt.
    # useful_words_count <= 130 or useful_words_count >= 70

    useful_words_count = 130

    if characters_count <= useful_words_count:

      confidence_sum = 0

      for confidence_value in res:

        # considering those words or characters whose accuracy in res array >= 98 %
        if confidence_value[2] >= 0.98:
          
          # Multiplying by thousand to increase magnitude of confidence_sum variable
          confidence_sum = confidence_sum + (confidence_value[2]*10000)

        else:

          confidence_sum = confidence_sum + confidence_value[2]
      # in this array I appended confidence_sum,character count values of Receipt at different angles
      confidence_and_angle_array.append((confidence_sum,angle,characters_count))
  
  #print(confidence_and_angle_array)
  return confidence_and_angle_array


def crop_receipt(img):
  """This function take image as input and remove skewness from Image
     And cropping image using 2 Algorithm approach inside function
     
     It is cropping Image using words detected by OCR and then fetching the outermost coordinates
     from all four corners of receipt.
     Using Histogram approach to deskew receipt image more accurately

  Parameters
  ---------

  img: array
    Image could be skewed or deskewed or can be of rotated at any angle

  Returns
  -------

  array:
    A 3D Array of Image,Image is cropped and deskewed

  """


  #Using easy OCR for reading english language and detecting coordinates of all the words
  # storing all detected words with bounding box coordinates and accuracy in array
  reader = easyocr.Reader(['en'])

  array_of_detected_characters = reader.readtext(img)


  # Initializing all four coordinates by reading coordinate values from zero index in result arary
  x_left = array_of_detected_characters[0][0][0][0]

  x_right = array_of_detected_characters[0][0][1][0]

  y_top =   array_of_detected_characters[0][0][0][1]

  y_bottom = array_of_detected_characters[0][0][3][1]

  # Iterating on each coordinate of characters detected by easy OCR
  # Saving coordinates of exterior most characters on all four sides(Imagine a Paper)
  for i in array_of_detected_characters:

    x_left = abs(min(x_left,i[0][0][0]))

    x_right = abs(max(x_right,i[0][1][0]))

    y_top = abs(min(y_top,i[0][0][1]))

    y_bottom = abs(max(y_bottom,i[0][3][1]))
    

  #print(int(x_left),(x_right),"-",int(y_top),int(y_bottom))
  cropped_image = img[int(y_top):int(y_bottom),int(x_left):int(x_right)] #Crop_Image

  #Passing cropped image to another fucntion to remove skewness in Image
  ans_arr = deskew_croped_receipt(cropped_image) # Crop_Image


  
  #  If our ans_arr can't detect any words on basis of hyperparameter
  #  I am using content driven image cropping inside if statement for more accuracy
  if len(ans_arr) == 0:

    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    denoise = cv2.fastNlMeansDenoising(gray,None,10,7,21)

    threshold_image = cv2.adaptiveThreshold(denoise,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,31,2)

    image_matrix = np.array(threshold_image)
    
    #inverting bgr value of white and black color
    # black == 1
    # white == 0
    inverse_image_pixel = 255-image_matrix

    # summing all the values of pixel count along y axis of an image
    y_axis_pixel_sum = inverse_image_pixel.sum(axis = 1)

    # Iterating from top to bottom along y axis of an image and storing coordinate of fist_pixel_occurence
    y_top_coordinate = 0
    for index, pixel in enumerate(y_axis_pixel_sum):
    
      if pixel >= 10:

        y_top_coordinate = index
        break

    # Iterating from bottom to top along y axis of an image and storing coordinate of fist_pixel_occurence
    y_bottom_coordinate = 0
    for ind, pix in reversed(list(enumerate(y_axis_pixel_sum))):
  
      if pix >= 10:

        y_bottom_coordinate = ind
        break

    # summing all the values of pixel count along x axis of an image
    x_axis_pixel_sum = inverse_image_pixel.sum(axis = 0)


    # Iterating from left to right along x axis of an image and storing coordinate of fist_pixel_occurence
    x_left_coordinate = 0
    for index, pixel in enumerate(x_axis_pixel_sum):

      if pixel >= 10:

        x_left_coordinate = index
        break


    # Iterating from right to left along x axis of an image and storing coordinate of fist_pixel_occurence
    x_right_coodinate= 0
    for ind, pix in reversed(list(enumerate(x_axis_pixel_sum))):

      if pix >= 10:

        x_right_coodinate = ind
        break
    

    # checking conditions to look for default skewed case(horizontal = 180, vertical = 90)
    if x_left_coordinate < y_top_coordinate:

      # rotating clockwise
      skew_change_default_90 = imutils.rotate_bound(cropped_image,90)

      # rotating anticlockwise
      skew_change_default_opposite_90 = imutils.rotate_bound(cropped_image,-90)

      #ouput array for temp_skew_optimization function
      best_array_output = []

      def temp_skew_optimization(image,angle):
        """Get image and rotated angle and give best accuracy of easy OCR on image with input rotated angle

        Parameters
        ----------

        img: array
           A 2D image array for detecting skewness
        angle: int
           An angle associated with image,image is default rotated at input angle

        Returns
        -------
           list
              A list of confidence( accuracy ) value at different angle of easy OCR
        """ 
        
        reader = easyocr.Reader(['en'])

        res = reader.readtext(image)

        confidence_sum_default = 0

        for confidence_value in res:

          if confidence_value[2] >= 0.98:

            confidence_sum_default = confidence_sum_default + (confidence_value[2]*10000)

          else:

            confidence_sum_default = confidence_sum_default + (confidence_value[2] * (-1000))

        best_array_output.append(confidence_sum_default)

        

      # calling temp_skew_optimization function to compare best rotation angle(90,-90)
      temp_skew_optimization(skew_change_default_90,90)
      temp_skew_optimization(skew_change_default_opposite_90,-90)

      # taking best index and angle
      best_rotation = best_array_output.index(max(best_array_output))
      #print(best_rotation)

      # taking best_rotation and returning best_rotated_image
      if best_rotation == 0:
        return skew_change_default_90
      elif best_rotation == 1:
        return skew_change_default_opposite_90
    
    # if every condition above not satisfy then return default cropped_image
    else:
      return cropped_image


  else:

    index = ans_arr.index(max(ans_arr))
    word_count_at_max_accuracy = ans_arr[index][2]

    if ans_arr[index][1] == 90:
      rotated_img = imutils.rotate_bound(cropped_image,90)
      return rotated_img

    # here checking if word_count at maximum accuracy is greater than default 0,90.270 degree case
    best_angle = 0
    flag = False
    for i in range(len(ans_arr)):
      if ans_arr[i][1] == 0 or ans_arr[i][1] == 90 or ans_arr[i][1] == 270 :
          if ans_arr[i][2] < word_count_at_max_accuracy:
            word_count_at_max_accuracy = ans_arr[i][2]
            best_angle = ans_arr[i][1]
            print("best angle is: ",best_angle)
            
            flag = True

    if flag:
      rotated_image = imutils.rotate_bound(cropped_image,best_angle)
      return rotated_image
    
    angle = max(ans_arr)[1]
    
    #Final deskewed and cropped image returning here
    rotated_img = imutils.rotate_bound(cropped_image,angle)
    return rotated_img

def preprocess(image_path):

  """This function call crop_receipt function
  for deskewness and cropping of manual clicked receipt image
  And preprocessing( converting image in binary form) image 

  Parameters
  ----------

  image_path: str
    The file location of image

  Returns
  -------
  Array
     A 2D array of deskew and cropped image converted into binary form
  """
  try:
    #Reading image using opencv and this convert image_path in array
    img = cv2.imread(image_path)
    #print(img)

    #calling crop_receipt function 
    crop_img = crop_receipt(img)

    #thresholding cropped Image
    gray = cv2.cvtColor(crop_img, cv2.COLOR_BGR2GRAY)
    denoise = cv2.fastNlMeansDenoising(gray, None, 10, 7, 21)
    threshold_image = cv2.adaptiveThreshold(denoise, 255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY, 31, 2)

    return threshold_image

  except:
    print("Image file extension should be suppported by CV2")


#main
image_path = "/content/sample_data/00EB3E618E9D433CA591-1.jpg"
cv2_imshow(cv2.imread(image_path))
deskewed_cropped_image = preprocess(image_path)
cv2_imshow(deskewed_cropped_image)
