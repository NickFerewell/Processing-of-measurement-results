from scipy import stats

print(stats.t.ppf(q=0.95, df=10))
print(stats.t.interval(confidence=0.95, df=9))