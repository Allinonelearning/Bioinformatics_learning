# 列表
# letters = ["A","b","C"]
# letters_uppercase = [letter.upper() for letter in letters]
# print(letters_uppercase)
# z = [0] * 100
# number = list(range(1,20,2))
# print(number)


# try:
#     with open("Data_Type.py") as file:
#         print("File opened.")
#     age =  int(input("Age:"))
#     xfactor = 10/age
# except (ValueError,ZeroDivisionError):
#     print("You didn't enter a valid age.")
# else:
#     print("No exceptions were thrown.")

# from timeit import timeit
# code1 = """
# def calculate_xfactor(age):
#     if age <= 0:
#         raise ValueError("Age cannot be 0 or less.")
#     return 10/age

# try:
#     calculate_xfactor(-1)
# except ValueError as error:
#     pass

# """
# print("first",timeit(code1, number=10000))


number = [1,2]
number.append(3)










