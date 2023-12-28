import random
import sys

if len(sys.argv) > 1:
    N = int(sys.argv[1])

    ans = [0 for _ in range(N)]
    b = [0 for _ in range(N)]
    table2 = [[0 for _ in range(N)] for _ in range(N)]

    for i in range(N):
        table2[i][i] = -2
        if i > 0:
            table2[i][i - 1] = 1
        if i < N - 1:
            table2[i][i + 1] = 1

    for i in range(N):
        ans[i] = random.random()  # Assigning random values to 'ans' between 0 and 1

    with open('matrix', 'w') as f:
        for i in range(N):
            for j in range(N):
                f.write(str(table2[i][j]))
                f.write('\t')
                b[i] += table2[i][j] * ans[j]
            f.write('\n')

    with open('matrix', 'a') as f:
        for i in range(N):
            f.write(str(b[i]))
            f.write('\n')

    with open('ans', 'w') as f:
        for i in range(N):
            f.write(str(ans[i]))
            f.write('\n')

else:
    print("Please input a number for matrix size")
