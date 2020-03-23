###### Simulated illness spread of COVID-19 in Sweden #####
###### (c) 2020 Jasmine Gardner, PhD ######



import time
import numpy as np 
from numpy import random 
import math
import pandas as pd
from matplotlib import pyplot as plt

f=open("out_file.txt", "w")
# Population
pop=1000

####### Initialize pandas array #############
pop_df=np.arange(0, pop, dtype=int)
# Initiate dataframe to hold people
population=pd.DataFrame(data=pop_df, columns=["Person"], dtype=int)
print((population))

###### Initialize ages ############
# Age ranges and distribution
# From SCB, population per month by region, age, and sex, National data in 5 year age intervals, Jan 2020, http://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__BE__BE0101__BE0101A/BefolkManad/

age_start=[0, 10, 20, 30, 40, 50, 60, 70, 80]
age_dist=[0.10721, 0.11553, 0.12439, 0.13475, 0.12578, 0.12704, 0.10782, 0.09952, 0.05796]

#Initialize age
age=[]
for i in range(0,pop):
	# Random age based on distribution
	age_ind=random.choice(np.arange(0, len(age_start)), p=age_dist)
	age.append(age_start[age_ind]+random.rand()*10)
#Merge with population data
population['Age']=age

###### Initialize Households #########
### According to OCED, Sweden has an average household size of approx 2 people. We can include a more representative distribution here if desired. 
### Each household needs an adult as a head of household.  Will first ignore children and place head of households, after each household has 1 person, I will go back and add in rest of people randomly.

household=np.ones((pop,), dtype=int)
household=household*-1
num_households=int(pop/2)
print(num_households)
HH=0
HH1=0

### Fill up households with a single adult first.
while HH<num_households:
	while population.iloc[HH1]['Age']<20:
		HH1+=1
	household[HH1]=HH
	HH1+=1
	HH+=1

### Add remaining people to a household.
for i in range(0, pop):
	if household[i]==-1:
		household[i]=math.floor(random.rand()*(num_households))

household_size=np.histogram(household, bins=np.arange(num_households))
household_n=household_size[0]
print(household_n)
household_size_dist=np.histogram(household_size[0], bins=np.arange(10))

# Add household to pandas array
population['Household']=household
print(population.dtypes)
print("household", time.time())


############## Designate region #############

# Assign households to region. May not be completely accurate with people per household but should be close with many households. 
# Counties: population data from SCB quarter 4 2019.
tot_pop=10327589
county=np.array(["Stockholm", "Uppsala", "Sodermanland", "Ostergotland", "Jonkoping", "Kronoberg", "Kalmar", "Gotland", "Blekinge", "Skane", "Halland", "Vasta Gotaland", "Varmland", "Orebro", "Vastmanland", "Dalarna", "Gavleborg", "Vasternorrland", "Jamtland", "Vasterbotten", "Norrbotten"]) 
pop_county=np.array([2377081, 383713, 297540, 465495, 363599, 201469, 245446, 59686, 159606, 1377827, 333848, 1725881, 282414, 304805, 275845, 287966, 287382, 245347, 130810, 271736, 250093])

# Confirm population per county matches total population.
tot_pop1=0
for i in pop_county:
	tot_pop1=tot_pop1+i
if (tot_pop!=tot_pop1):
	print("POPULATIONS DO NOT MATCH", tot_pop, tot_pop1)

# Randomly assign household based on population distribution.
pop_dist=pop_county/tot_pop

HH_county=[]
for i in range(0, num_households):
	HH_county.append(random.choice(np.arange(0, len(pop_dist)), p=pop_dist))
county_pop=np.zeros(pop, dtype=int)

for i in range(0, pop):
	county_pop[i]=HH_county[int(population.iloc[i]['Household'])]

population['County']=county_pop

#Check population distribution
county_size=np.histogram(county_pop, bins=np.arange(len(county)), density=True)


######### Initialize job/school situation #############
#  Job types are 1 for preschool (age 1-5) (current participation in Sweden is 86% for ages 1-3 and 95% for ages 3-5 so we will assume 90% participate in a preschool program), 2 for elementary school (age 6-15) (assume 100% participation), 3 for high school/university (age 15-22) (assume 100% participation), and 4 for a job (age 22-75) of which 73.4% participate in workforce (from SCB).  Infants under 1 and elderly people over 75 do not go to a job.


job_status=np.zeros(pop, dtype=int)
job_status_county=np.zeros((5, len(county)), dtype=int) 
for i in range(0, pop):
	if (age[i]<1 or age[i]>76):
		job_status[i]=0
		job_status_county[0][int(county_pop[i])]+=1
	elif (age[i]>=1 and age[i]<6):
		if (random.rand()<0.900):
			job_status[i]=1
			job_status_county[1][int(county_pop[i])]+=1
		else:
			job_status[i]=0
			job_status_county[0][int(county_pop[i])]+=1
	elif (age[i]>=6 and age[i]<15):
		job_status[i]=2
		job_status_county[2][int(county_pop[i])]+=1
	elif (age[i]>=15 and age[i]<22):
		job_status[i]=3
		job_status_county[3][int(county_pop[i])]+=1
	elif (age[i]>=22 and age[i]<76):
		if (random.rand()<0.734):
			job_status[i]=4
			job_status_county[4][int(county_pop[i])]+=1
		else:
			job_status[i]=0
			job_status_county[0][int(county_pop[i])]+=1

#Check population distribution
job_size=np.histogram(job_status, bins=np.arange(5), density=True)

population['Job_Status']=job_status


########## Initialize workplace/school ############

#Assume (ppschool) people per school and (ppwork) per close workgroup. Can set up array and do mapped function when we need to specify different group sizes per type of job.
#OCED states 19.2 pupils per elementary school class and 21.2 in secondary class in Sweden. Assuming 20 pupils per class and 6 classes per school, average school size will be 120 students.
#"Workplace" average size will be 15 persons.  The actual distribution for micro (0-9), small(10-49), medium(50-249), and large(249+) is [0.9455, 0.04510, 0.0079, 0.0014].  We will assume that the average person comes in close contact with 15 coworkers during their normal workday; this could indicate the smaller group of core contacts in a larger business or outside contacts including regular deliveries and other support staff in micro-businesses. 
ppwork=10
ppschool=120
num_workplaces=np.ceil(job_status_county/int(ppschool))
num_workplaces[4]=np.ceil(job_status_county[4]/int(ppwork))
max_num_workplaces=np.amax(num_workplaces)

#Workplace array
workplace=np.zeros(pop, dtype=int)
workplace_num=np.zeros((5, len(county), int(max_num_workplaces))) 
for i in range(0, len(county)):
	for index, row in population[population.County==i].iterrows():
		workplace[index]=np.floor(random.rand()*num_workplaces[int(job_status[index])][int(county_pop[index])])
		workplace_num[int(job_status[index])][int(county_pop[index])][int(workplace[index])]+=1;
population['Workplace']=workplace


####### Seed initial infection #######
## Infections by county
initial_infections=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
symptomatic=0.33 # Percent which are believed to be symptomatic.  Those with symptoms are assumed to be twice as infectious as asymptomatic carriers. (Ferguson Nature 2020)

infections=np.zeros(pop)
severe=np.zeros(pop) #50% of cases are severe.  Determines infectiousness. 
tau=np.ones(pop)*-100 #Onset of infection
sympt=np.zeros(pop) 

for i in range(0, len(initial_infections)):
	infected_init=0
	person_infected=0
	while (infected_init<initial_infections[i] and person_infected<pop):
		if (population.iloc[person_infected]['County']==i):
			infections[person_infected]=1
			severe[person_infected]=round(random.rand())
			if (random.rand()<symptomatic):
				sympt[person_infected]=1
			tau[person_infected]=0
			infected_init+=1
			
		person_infected+=1
		

population['Infected']=infections

print(time.time())

#### Start simulation ####


##### Initialize constants ######
betah=0.627 # Scaled from betah=0.4 in influenza pandemic with R0=1.6, COVID-19 R0=2.4 (Ferguson 2020)
betac=0.1 # Scaled from betah=0.075 in influenza pandemic with R0=1.6, COVID-19 R0=2.4 (Ferguson 2020)
betap=[0.0, 1.254, 1.254, 1.254, 0.627] # Spread in all types of schools (preschool to college) is twice that of 
kappa=1 #What is this parameter??
alpha=0.8 # From Ferguson Nature 2006
omega=2 # From Ferguson Nature 2006
psi=np.array([0.0, 0.1, 0.2, 0.25, 0.5]) # Accounts for absenteeism based on severe illness. Ferguson Nature 2006
zeta=np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.50, 0.25, 0.25, 0.25])   # Travel related parameter for community transmission. Ferguson Nature 2006
#Leaving out rho from Ferguson 2006.  This is a measure of how infectious person is.  For now, we will assume all people are the same.

### Rates of hospitalization and death by age group as specified by age start as percent of symptomatic cases.
hosp_symptomatic=np.array([0.001, 0.003, 0.012, 0.032, 0.049, 0.102, 0.166, 0.243, 0.273]) # From Ferguson 2020
icu=np.array([0.05, 0.05, 0.05, 0.05, 0.063, 0.122, 0.274, 0.432, 0.709]) # percent of hospitalized that need icu from Ferguson 2020.
fatal_in_icu=0.5
fatal_symptomatic=np.array([0.00002, 0.00006, 0.0003, 0.0008, 0.0015, 0.006, 0.022, 0.051, 0.093])


# Initialize hospitalization and death for use to minimize calculation.
hosp=hosp_symptomatic*symptomatic
fatal=fatal_symptomatic*symptomatic
t=0
infect=np.zeros(pop)
recovered=np.zeros(pop)
hosp_pop=np.zeros(pop)
icu_pop=np.zeros(pop)
death=np.zeros(pop)
recovered=np.zeros(pop)

population['Hospitalized']=hosp_pop
population['ICU']=icu_pop
population['Tau']=tau
population['Death']=death
population['Recovered']=recovered
population['Symptomatic']=sympt

print("Currently Infected", "Death", "HOSP", "ICU", "Recovered", "Susceptible", file=f, flush=True)
for t in range(5, 2000):
	start=time.time()
	print("Timestep ", t)
#	print("Timestep ", t, population[(population.Infected == 1) & (population.Recovered == 0)])
	###############First step, find infectivity ######
	
	#### Only Susceptible people can get the virus.
	for index1, row1 in population[population.Infected==0].iterrows():
		
		#### Only infected people can spread the virus.  Only include those that are currently infectious and not hospitalized.
		for index2, row2 in population[(population.Tau>(t-11.1)) & (population.Tau<(t-4.6)) & (population.Hospitalized==0) & (population.ICU==0)].iterrows():

			#### NOTE: Infectiousness is considered to be equal for all persons with a value of 1.

			#### Only consider people in the same county.
			if (county_pop[index1]==county_pop[index2]):
				###Determine kappa for infected person.  This is the infectiousness of the person based on time since infection started.  Latency period is 4.6 days.  Infection starts at 5.1 days and lasts for 6 days.  Sympotmatic people are twice as likely to infect others as asymptomatic.
				### NOTE: In original study, this a function based on time.  This can be added later. 
				if (t-tau[index2]<4.6):
					kappa=0
				elif (t-tau[index2]>11.1):
					kappa=0 # Recovered or dead
				elif (sympt[index2]==1):
					kappa=1
				else:
					kappa=0.5
				
				##### Household transmission #####
				if (household[index1]==household[index2]):

		#			print("In home")	
					household_tmp=int(household[index1])
					#### Account for household transmission. #######
					infect[index1]+=1*betah*kappa*(1+severe[index2]*(omega-1))/(pow(household_n[household_tmp-1],alpha));
					
				##### Workplace/School transmission ######
				### People must be in same workplace and job type. ####
				if (workplace[index1]==workplace[index2]) and (job_status[index1]==job_status[index2]) :
		#			print("In workplace")
					#### Account for workplace transmission. #######
					infect[index1]+=(1*betap[int(job_status[index1])]*kappa*(1+severe[index2]*(omega*psi[int(job_status[index1])]-1))/workplace_num[int(job_status[index1])][int(county_pop[index1])][int(workplace[index1])]);

				###### Community transmission #####
				# Replaced distance metric with overall population similar to workplace transmission. 

				age_group=math.floor(age[index1]/5)
				infect[index1]+=1*zeta[age_group]*betac*kappa*(1+severe[index2]*(omega-1))/(pop_county[int(county_pop[index1])]);


		### Probability of being infected ####
		if (infect[index1]>0):
			infect_prob=(1-math.exp(-infect[index1]))

			### Monte carlo type prediction ###
			if (random.rand()<infect_prob):
	#			print("infected person", infect_prob)
				infections[index1]=1
				tau[index1]=t
				severe[index1]=round(random.rand())
				sympt[index1]=int(random.rand()<symptomatic)

	
	#### After 5 days of symptoms (10.1 days after infection), randomly place people in hospital based on age.  Also determine those who need ICU care. Note: Only symptomatic people can go to the hospital.  ######
	for index, row in population[(population['Tau']==(t-10)) & (population['Symptomatic']==1)].iterrows():
		age_group=math.floor(age[index]/10)

		if (random.rand()<hosp[age_group]*symptomatic):
			if (random.rand()<icu[age_group]):
				icu_pop[index]=1	
			else:
				hosp_pop[index]=1

	#### 15 days after onset of symptoms, people randomly die based on age based percentages.  ICU patients are removed from ICU into regular hospital. ####
	for index, row in population[(population['Tau']==(t-20)) & (population['Symptomatic']==1)].iterrows():
		age_group=math.floor(age[index]/10)
		if (icu_pop[index]==1):
			icu_pop[index]=0
			death[index]=round(random.rand()) # 50% chance of death in ICU
			if (death[index]==0):
				hosp_pop[index]=2  # Move from ICU to regular hospital.
		elif (random.rand()<fatal[age_group]):
			death[index]=1
			hosp_pop[index]=0

	#### Release from hospital   ######	
	for index, row in population[(population['Tau']>=(t-18)) & (population['Hospitalized']==1)].iterrows():
		#### After 8 days in hospital, hospitalized patients are released from hospital. ####
		if (tau[index]==t-18):
			hosp_pop[index]=0
			recovered[index]=1
	for index, row in population[(population['Tau']>=(t-26)) & (population['Hospitalized']==2)].iterrows():
		#### After 6 days in hospital after ICU stay, ICU patients released from hospital. ####
		if (tau[index]==t-26):
			print("here")
			hosp_pop[index]=0
			recovered[index]=1

	#### Recovered after 11 days (6 days of symptoms) if not in hospital/ICU. ####
	for index, row in population[(population['Tau']==(t-11)) & (population['ICU']==0) & (population['Hospitalized']==0) & (population['Infected']==1) & (population['Recovered']==0)].iterrows():
		recovered[index]=1



	population.Infected=infections
	population.Tau=tau
	population.Hospitalized=hosp_pop
	population.ICU=icu_pop
	population.Death=death
	population.Recovered=recovered
	population.Symptomatic=sympt

#	print(population[(population.Tau<(t-26)) & (population.Recovered==0) & (population.Tau>-1)])
#	print(population[(population.Death==1) & (population.Infected==0) & (population.Tau>-1)])
	end=time.time()
	print("time", end-start)
	print(t, len(population[(population.Tau>=(t-11)) & (population.Infected==1)]), len(population[population.Death==1]), len(population[population.Hospitalized==1]), len(population[population.ICU==1]), len(population[population.Recovered==1]), len(population[population.Infected==0]), file=f, flush=True)	
