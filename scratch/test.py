
print(lambda x: x + 1)
print((lambda x: x + 1)(2))

# use lambda to define a function: if true return x + 1, else return x - 1
print((lambda x: x + 1 if x > 0 else x - 1)(True))

lb = lambda x: x + 1 if x > 0 else x - 1
print(lb(True))

import logomaker as lm
