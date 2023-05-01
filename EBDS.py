import numpy as np
import matplotlib.pyplot as plt
import random

num_sims = 10
mut = 0.05
dose = 0#.6

g = 0.02
phi = 0.02
K = 100
lam = 1
b = 2
r = 0.3
s = 0.6

tend = 1000

sens_macro = []
resist_macro = []
strat_macro = []
t_macro = []
extinct = []

def death_drug(current_strat,m1): #fitness analogous to death due to drug in our case
        return m1/(lam+b*current_strat)

def run_this():
        sens = [10]
        resist = [0]
        strat = [0]
        t = [0]

        while t[-1] < tend:

                if t[-1] > 200 and t[-1] < 800:
                        m1 = dose
                else:
                        m1 = 0

                current_sens = sens[-1]
                current_resist = resist[-1]
                current_strat = strat[-1]

                #broken into birth, death, switching

                rates = [current_sens*s, current_sens*(m1/(lam+b*current_strat)                     +s*((current_sens+current_resist)/K)), current_sens*g, current_resist*r, current_resist*(current_sens+current_resist)/K, current_resist*phi]

                rate_sum = sum(rates)

                if rate_sum == 0:
                        extinct.append(t[-1])
                        break

                tau = np.random.exponential(scale=1/rate_sum)

                t.append(t[-1] + tau)

                rand = random.uniform(0,1)

                #Sensitive cell division event
                if rand * rate_sum <= rates[0]:
                        resist.append(resist[-1])
                        if random.uniform(0,1)>mut:
                                strat.append(strat[-1])
                                sens.append(sens[-1] + 1)
                        else:
                                new_strat = strat[-1]+np.random.normal(0,0.01)
                                if death_drug(new_strat,m1)<death_drug(strat[-1],m1):
                                        strat.append(new_strat)
                                        sens.append(sens[-1] + 1)
                                elif death_drug(new_strat,m1)==death_drug(strat[-1],m1):
                                        strat.append(strat[-1])
                                        sens.append(sens[-1] + 1)
                                else:
                                        strat.append(strat[-1])
                                        sens.append(sens[-1])

                #Sensitive cell death event
                elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):
                        sens.append(sens[-1] - 1)
                        resist.append(resist[-1])
                        strat.append(strat[-1])

                #Sensitive switch event
                elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):
                        sens.append(sens[-1]-1)
                        resist.append(resist[-1] + 1)
                        strat.append(strat[-1])

                #Resistant division event
                elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):
                        sens.append(sens[-1])
                        resist.append(resist[-1] + 1)
                        strat.append(strat[-1])

                #Resistant death event
                elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):
                        sens.append(sens[-1])
                        resist.append(resist[-1]-1)
                        strat.append(strat[-1])

                #Resistant switch event
                elif rand * rate_sum > sum(rates[:5]) and rand * rate_sum <= sum(rates[:6]):
                        resist.append(resist[-1]-1)
                        sens.append(sens[-1] + 1)
                        strat.append(strat[-1])

        sens_macro.append(sens)
        resist_macro.append(resist)
        strat_macro.append(strat)
        t_macro.append(t)

for i in range(num_sims):
        run_this()

lwi = 0.3

extinct_events = len(extinct)
avg_extinct = np.mean(extinct)
std_extinct = np.std(extinct)

print(extinct_events)
print(avg_extinct)
print(std_extinct)

plt.figure()
plt.subplot(211)
plt.title('Control')
plt.plot(t_macro[0],sens_macro[0], label='Sensitive', c='r',lw=lwi)
plt.plot(t_macro[0],resist_macro[0], label='Resistant', c='b',lw=lwi)
for i in range(1,num_sims):
    plt.plot(t_macro[i],sens_macro[i], c='r',lw=lwi)
    plt.plot(t_macro[i],resist_macro[i], c='b',lw=lwi)
plt.legend()
plt.grid(True)
ax = plt.gca()
#ax.axvspan(200, 800, facecolor=(0.5,0.5,0.5))
plt.ylabel('Population Size')
plt.subplot(212)
plt.plot(t_macro[0],strat_macro[0], c='k',lw=lwi)
for i in range(1, num_sims):
    plt.plot(t_macro[i],strat_macro[i], c='k',lw=lwi)
plt.grid(True)
ax = plt.gca()
#ax.axvspan(200, 800, facecolor=(0.5,0.5,0.5))
plt.xlabel('Time')
plt.ylabel('Drug Resistance')
plt.tight_layout()
plt.show()