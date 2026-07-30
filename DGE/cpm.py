import csv

att1_sum = 45809623
att2_sum = 42528063
att3_sum = 49742697
swim1_sum = 34300954
swim2_sum = 47969916
swim3_sum = 42620690

with open("counts_dup_prim_7-28.csv", mode='r', newline='') as f, open("norm_counts_7-28.csv", mode='w') as n:
    reader = csv.reader(f)
    writer = csv.writer(n)
    header = next(reader)
    writer.writerow(header)
    for row in reader:
        cpm_att1 = (float(row[1])/att1_sum) * 1000000
        cpm_att2 = (float(row[2])/att2_sum) * 1000000
        cpm_att3 = (float(row[3])/att3_sum) * 1000000
        cpm_swim1 = (float(row[4])/swim1_sum) * 1000000
        cpm_swim2 = (float(row[5])/swim2_sum) * 1000000
        cpm_swim3 = (float(row[6])/swim3_sum) * 1000000
        new_row = [row[0], cpm_att1, cpm_att2, cpm_att3, cpm_swim1, cpm_swim2, cpm_swim3]
        writer.writerow(new_row)



