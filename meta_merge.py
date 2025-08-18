import pandas as pd
import numpy as np
import sys

# 取得 pheno 參數
pheno = sys.argv[1] if len(sys.argv) > 1 else "pheno"

# 讀入 fixed 與 random 結果
fixed = pd.read_csv('meta_fixed_results.txt', sep='\t')
random = pd.read_csv('meta_random_results.txt', sep='\t')

# 合併
merged = pd.merge(
    fixed[['ID', 'META_LOG_OR', 'META_SE', 'META_P', 'Q', 'I2']],
    random[['ID', 'META_LOG_OR', 'META_SE', 'META_P', 'Q', 'I2']],
    on='ID',
    suffixes=('_fixed', '_random')
)

# 新增 sensitivity analysis 欄
def sensitivity(row):
    pf, pr = row['META_P_fixed'], row['META_P_random']
    i2 = row['I2_fixed']
    logor_f = row['META_LOG_OR_fixed']
    logor_r = row['META_LOG_OR_random']

    if i2 < 30:
        return 'Use fixed-effect'
    elif 30 <= i2 < 50:
        return 'Moderate heterogeneity'
    elif i2 >= 50:
        if (pf < 0.05 and pr > 0.05) or (pf > 0.05 and pr < 0.05) or (np.sign(logor_f) != np.sign(logor_r)):
            return 'Use random-effect (conflicting results)'
        else:
            return 'High heterogeneity but consistent results'
    else:
        return 'Check input'

merged['sensitivity_analysis'] = merged.apply(sensitivity, axis=1)

# 輸出
output_file = f"{pheno}_meta_merge_results.txt"
merged.to_csv(output_file, sep='\t', index=False)
print(output_file)