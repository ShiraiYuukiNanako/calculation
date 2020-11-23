import csv
import mat

hole = list()
NS, EW, GRADE = dict(), dict(), dict()

with open('kriging.csv', newline='') as csvfile:
    file = csv.reader(csvfile)
    for row in file:
        hole.append(row[0])
        EW[row[0]] = row[1]
        NS[row[0]] = row[2]
        GRADE[row[0]] = row[3]

NS.pop('HOLE')
EW.pop('HOLE')
GRADE.pop('HOLE')
hole.remove('HOLE')


def NS_count_pair(lag_distance):
    result = []
    for i in range(len(hole)):
        for j in range(len(hole)):
            if int(NS[hole[i]]) - int(NS[hole[j]]) == lag_distance and int(EW[hole[i]]) == int(EW[hole[j]]):
                result.append((hole[i],hole[j]))
    return result


def EW_count_pair(lag_distance):
    result = []
    for i in range(len(hole)):
        for j in range(len(hole)):
            if int(EW[hole[i]]) - int(EW[hole[j]]) == lag_distance and int(NS[hole[i]]) == int(NS[hole[j]]) :
                result.append((hole[i],hole[j]))
    return result


lag_distance = [200, 400, 600, 800, 1000]

NS_count, EW_count = dict(), dict()

for lag in lag_distance:
    NS_count[lag] = len(NS_count_pair(lag))
    EW_count[lag] = len(EW_count_pair(lag))

print(NS_count)
print(EW_count)


def gamma(h):
    grade_sum = 0
    for pair in EW_count_pair(h):
        grade_sum += mat.pow(float(GRADE[pair[0]]) - float(GRADE[pair[1]]), 2)
    return grade_sum/(2*len(NS_count_pair(h)))


gammas = dict()

for lag in lag_distance:
    gammas[lag] = gamma(lag)
print(gammas)