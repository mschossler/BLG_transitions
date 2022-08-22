import csv

f_Up = open('densitiesUp.csv', mode='r')
csv_reader_Up = csv.reader(f_Up, delimiter=',')

f_Um = open('densitiesUm.csv', mode='r')
csv_reader_Um = csv.reader(f_Um, delimiter=',')

filling_order_Upositive = []
for row in csv_reader_Up:
    next_state_to_be_filled = [el for el in row if el not in filling_order_Upositive]
    filling_order_Upositive.extend(next_state_to_be_filled)

filling_order_Unegative = []
for row in csv_reader_Um:
    next_state_to_be_filled = [el for el in row if el not in filling_order_Unegative]
    filling_order_Unegative.extend(next_state_to_be_filled)

old_filling_order_Upositive = ['LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sdown',
                               'LL1_Kp_Sdown', 'LL1_Km_Sup', 'LL1_Kp_Sup', 'LL2_Kp_Sdown', 'LL2_Kp_Sup']

old_filling_order_Unegative = ['LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LL0_Kp_Sdown', 'LL0_Km_Sdown', 'LL0_Kp_Sup', 'LL0_Km_Sup', 'LL1_Kp_Sdown',
                               'LL1_Km_Sdown', 'LL1_Kp_Sup', 'LL1_Km_Sup', 'LL2_Km_Sdown', 'LL2_Km_Sup']

print(filling_order_Unegative)
print(old_filling_order_Unegative)

print(old_filling_order_Upositive == filling_order_Upositive)
print(old_filling_order_Unegative == filling_order_Unegative)
