
# Fashion MNIST Classifier using CNN and Pre-trained ResNet50

## Overview
This repository contains two implementations for classifying the Fashion MNIST dataset:

1. **Basic CNN Model**: A basic Convolutional Neural Network (CNN) model for image classification.
2. **Pre-trained ResNet50 Model**: A more advanced model that utilizes the pre-trained ResNet50 architecture, fine-tuned for the Fashion MNIST dataset.

The goal of this project is to improve the validation accuracy of the basic CNN model to above 90% by using a more sophisticated architecture, namely ResNet50.

## Dataset
The dataset used in this project is the Fashion MNIST dataset, which consists of 60,000 grayscale images of fashion items, each of size 28x28 pixels. There are 10 classes in total, including T-shirts, trousers, pullovers, dresses, and more.

## Models

### 1. Basic CNN Model
The first model implemented is a simple multi-layer perceptron-like CNN architecture that includes:
- Flattening of the input images from 2D (28x28) to 1D.
- Two dense layers:
  - Dense Layer 1: ReLU activation.
  - Dense Layer 2: Softmax activation for the output layer.

This model achieves around 87% accuracy.

```python
# Basic CNN Model Example
model = Sequential()
model.add(Flatten(input_shape=(img_width, img_height, 1)))
model.add(Dense(num_classes, activation="relu"))
model.add(Dense(num_classes, activation="softmax"))
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
```

### 2. Pre-trained ResNet50 Model
The second model is based on the ResNet50 architecture, which is a deep residual network that has been pre-trained on the ImageNet dataset. By leveraging transfer learning, we can fine-tune ResNet50 to perform well on the Fashion MNIST dataset, even though it was originally trained on a completely different dataset.

#### Key Modifications:
- The input images are resized from 28x28 to 224x224 to match the expected input size of ResNet50.
- The grayscale images are converted to RGB by repeating the single channel three times.
- The ResNet50 base model is used with pre-trained weights from ImageNet, and we freeze the base model to prevent its weights from being updated during training.
- We add a custom fully connected layer with a Dense(512) layer, followed by a dropout layer to reduce overfitting, and finally an output softmax layer for classification into 10 classes.

This model achieves a validation accuracy of above 90%, outperforming the basic CNN model.

## How to Run

### 1. Basic CNN Model
To run the basic CNN model, follow these steps:

1. Clone the repository and navigate to the project directory.
2. Install the required dependencies using:

```bash
pip install tensorflow wandb
```

3. Run the `cnn_basic.py` file:

```bash
python cnn_basic.py
```

### 2. Pre-trained ResNet50 Model
To run the pre-trained ResNet50 model, follow these steps:

1. Clone the repository and navigate to the project directory.
2. Install the required dependencies using:

```bash
pip install tensorflow wandb
```

3. Run the `resnet50_fashion.py` file:

```bash
python resnet50_fashion.py
```

## Results

| Model                  | Accuracy |
|------------------------|----------|
| **Basic CNN Model**     | ~87%     |
| **Pre-trained ResNet50**| >90%     |

By leveraging the pre-trained ResNet50 model, we achieved a significant improvement in the validation accuracy compared to the basic CNN model.

## Key Techniques Used
- **Transfer Learning**: We utilized the powerful ResNet50 architecture pre-trained on the ImageNet dataset, which helped improve classification performance on the Fashion MNIST dataset.
- **Data Preprocessing**: The grayscale images were reshaped and resized to fit the input requirements of the ResNet50 model.
- **Dropout Regularization**: A dropout layer was added to reduce overfitting in the final dense layers.

## Future Work
- **Fine-tuning the ResNet50 base layers**: Instead of freezing all ResNet layers, we could fine-tune the last few layers of the ResNet model for further improvements.
- **Experimenting with other architectures**: Testing more advanced pre-trained models like EfficientNet or InceptionV3 may yield even better results.
