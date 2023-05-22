import numpy as np
inputSize = 1024
level=int(np.log2(inputSize))

Wn_re = np.zeros(inputSize)
Wn_im = np.zeros(inputSize)

for i in range(inputSize):
    Wn_re[i] = np.cos(-2 * np.pi * i / inputSize)
    Wn_im[i] = np.sin(-2 * np.pi * i / inputSize)


Wn_new_re = np.zeros((level, inputSize//2))
Wn_new_im = np.zeros((level, inputSize//2))


m = inputSize
for i in range(int(level)):
    m = m // 2
    N = inputSize // m
    new_index = 0
    for p in range(m):
        for j in range(N//2):
            idx = j * m
            Wn_new_re[i][new_index] = Wn_re[idx]
            Wn_new_im[i][new_index] = Wn_im[idx]
            new_index += 1

with open("./Wn_new.txt", "w") as f:
    for i in range(level):
        f.write(f"\n\nec::Float Wn_new_re_{i}[WINDOW_SIZE/2] = " + "{")
        for j in range(inputSize//2):
            f.write(str(Wn_new_re[i][j]) + 'f')
            if j != (inputSize//2-1):
                f.write(", ")
            else:
                f.write("};")

    f.write("\n\n\n\n")
    f.write("*"*200)
    f.write("\n\n\n\n")

    for i in range(level):
        f.write(f"\n\nec::Float Wn_new_im_{i}[WINDOW_SIZE/2] = " + "{")
        for j in range(inputSize//2):
            f.write(str(Wn_new_im[i][j]) + 'f')
            if j != (inputSize//2-1):
                f.write(", ")
            else:
                f.write("};")

