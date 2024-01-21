import numpy as np

class LinearSVM:
    def __init__(self, learning_rate=0.01, num_epochs=1000, C=1.0):
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.C = C
        self.W = None
        self.b = None

    def fit(self, X, y):
        samples, features = X.shape
        y1 = np.where(y <= 0, -1, 1)

        # Initialize weights and bias
        self.W = np.zeros(features)
        self.b = 0

        # Gradient descent
        for epoch in range(self.num_epochs):
            for i in range(samples):
                if y1[i] * (np.dot(X[i], self.W) - self.b) >= 1:
                    
                    self.W -= self.learning_rate * (2 * self.C * self.W)
                else:
                    # Update the weights and bias
                    self.W -= self.learning_rate * (2 * self.C * self.W - np.dot(X[i], y1[i]))
                    self.b -= self.learning_rate * y1[i]

    def predict(self, X):
        return np.sign(np.dot(X, self.W) - self.b)


if __name__ == '__main__':


    np.random.seed(98)
    X = np.random.randn(300, 2)
    y = np.where(X[:, 0]**2 + X[:, 1]**2 > 1, 1, -1)

    # Create and train the SVM
    svm = LinearSVM(learning_rate=0.01, num_epochs=1000, C=1.0)
    svm.fit(X, y)

    # Predict and print accuracy
    y_pred = svm.predict(X)
    accuracy = np.mean(y_pred == y)
    print(f"Accuracy: {accuracy * 100:.2f}%")
