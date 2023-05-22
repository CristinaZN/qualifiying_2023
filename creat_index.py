import numpy as np

inputSize=8
level = int(np.log2(inputSize))
m = inputSize
new_index = []
# original_index = [i for i in range(inputSize)]
# last_even_groups = [original_index]
# ast_odd_groups = []
# for l in range(1,level, 1):
#     num_of_group = int(pow(2, l+1))
#     even_groups = []
#     odd_groups = []
#     for i in range(num_of_group):
#         if i % 2 == 0:
#             even_groups.append()
#             pass
#         if i % 2 == 1:
#             pass

m = inputSize
for l in range(level):
    x=[]
    y=[]
    m = m // 2
    for i in range(inputSize):
        if i not in y:
            x.append(i)
            y.append(i+m)
    new_index.append(x)




m = inputSize
# Wn_re = np.zeros((level, inputSize//2))
# Wn_im = np.zeros((level, inputSize//2))
# for l in range(level):
#     m = inputSize // pow(2, l+1)
#     for i in range(inputSize//2):
#         Wn_re[l][i] = np.cos(-2 * np.pi * ((i + m)//(2 * m)) / pow(2, l+1))
#         Wn_im[l][i] = np.sin(-2 * np.pi * ((i + m)//(2 * m)) / pow(2, l+1))

# validate





with open("./fft_index.txt","w") as f:
    for l in range(level):
        f.write(f"int new_index_{l}[WINDOW_SIZE/2] = " + "{")
        for j in range(inputSize//2):
            f.write(str(new_index[l][j]))
            if j!=(inputSize//2 - 1):
                f.write(", ")
            else:
                f.write("};")
        f.write("\n")
    f.write("std::vector<int*> new_index = {\n")
    for l in range(level-1, 0, -1):
        f.write(f"new_index_{l}")
        if l != 0:
            f.write(", ")
        else:
            f.write("\n};")
    f.write("\n")

# with open("./Wn_new.txt", "w") as f:
#     for l in range(level):
#         f.write(f"\n\nec::Float Wn_re_{l}[WINDOW_SIZE/2] = " + "{")
#         for i in range(inputSize//2):
#             f.write(str(Wn_re[l][i]) + 'f')
#             if i != (inputSize//2-1):
#                 f.write(", ")
#             else:
#                 f.write("};")

#     f.write("\n\n\n\n")
#     f.write("*"*200)
#     f.write("\n\n\n\n")

#     for l in range(level):
#         f.write(f"\n\nec::Float Wn_im_{l}[WINDOW_SIZE/2] = " + "{")
#         for i in range(inputSize//2):
#             f.write(str(Wn_im[l][i]) + 'f')
#             if i != (inputSize//2-1):
#                 f.write(", ")
#             else:
#                 f.write("};")