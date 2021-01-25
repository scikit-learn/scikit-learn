"""
Write a program that prints the integers from 1 to 100 (inclusive).
But:
for multiples of three, print Fizz (instead of the number)
for multiples of five, print Buzz (instead of the number)
for multiples of both three and five, print FizzBuzz (instead of the number)
"""

# Author: Danick Daloya

for i in range(1, 101):
    if i % 15 == 0:
        print("FizzBuzz")
    elif i % 3 == 0:
        print("Fizz")
    elif i % 5 == 0:
        print("Buzz")
    else:
        print(i)