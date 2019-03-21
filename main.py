#last updated 4/7/18
from sympy.solvers import solve
from sympy import Symbol, cos, sin, pi, sqrt, sympify
from decimal import Decimal
from termcolor import colored,cprint
from itertools import combinations
import replit
replit.clear()
#welcome message
cprint('Welcome to my Physics Equations Solver.\nIt can solve all Physics I & II equations.\nIt can also string together equations.','yellow')

#ALL THE VARIABLES (variable list)
vl=[
  #CONSTANTS
  {
    'c':'speed of light',
    'g':'acceleration due to gravity',
    'eo':'permittivity constant',
    'uo':'permeability constant',
    'ke':"coulumb's constant",
    'G':'universal gravitational constant',
    'hc':'planck constant',
    'kb':'boltzmann constant',
    'Rc':'gas constant'
  },
  #MECHANICS
  {
    'A':'amplitude',
    'a':'acceleration',
    'aa':'angular acceleration',
    'd':'displacement',
    'E':'energy',
    'F':'force',
    'f':'frequency',
    'I':'rotational inertia',
    'K':'kinetic energy',
    'k':'spring constant',
    'L':'angular momentum',
    'l':'length',
    'm':'mass',
    'm1':'mass1',
    'm2':'mass2',
    'P':'power',
    'p':'momentum',
    'r':'radius or seperation',
    'T':'period',
    't':'time',
    'tau':'torque',
    'theta':'angle',
    'U':'potential energy',
    'u':'coefficient of friction',
    'v':'velocity (initial)',
    'vf':'velocity final',
    'W':'work done on system',
    'w':'angular velocity (initial)',
    'wf':'angular velocity final',
    'x':'position',
    'y':'height'
  },
  #ELECTRICITY AND MAGNETISM
  {
    'A':'area',
    'B':'magnetic field',
    'C':'capacitance',
    'd':'distance',
    'E':'electric field',
    'F':'force',
    'I':'current',
    'l':'length',
    'P':'power',
    'Q':'charge',
    'q':'point charge',
    'q1':'point charge 1',
    'q2':'point charge 2',
    'R':'resistance',
    'r':'seperation',
    't':'time',
    'U':'potential energy(stored)',
    'V':'electric potential/emf',
    'v':'speed',
    'kappa':'dielectric constant',
    'p':'resistivity',
    'theta':'angle',
    'flux':'flux'
  },
  #FLUID MECHANICS AND THERMAL PHYSICS
  {
    'A':'area',
    'A1':'area1',
    'A2':'area2',
    'F':'force',
    'h':'height/depth',
    'h1':'height1',
    'h2':'height2',
    'k':'thermal conductivity',
    'K':'kinetic energy',
    'L':'thickness',
    'm':'mass',
    'n':'number of moles',
    'N':'number of molecules',
    'P':'pressure',
    'Pi':'pressure initial',
    'P1':'pressure1',
    'P2':'pressure2',
    'p':'density',
    'Q':'heat',
    'T':'temperature',
    't':'time',
    'U':'internal energy',
    'V':'volume',
    'v':'speed',
    'v1':'speed1',
    'v2':'speed2',
    'W':'work done on a system',
  },
  #MODERN PHYSICS, WAVES, AND OPTICS
  {
    'd':'separation/distance',
    'di':'image distance',
    'do':'object distance',
    'E':'energy',
    'f':'frequency/focal length',
    'h':'height',
    'hi':'image height',
    'ho':'object height',
    'K':'kinetic energy',
    'L':'distance',
    'wave':'wavelength',
    'M':'magnification',
    'm':'mass',
    'N':'an integer',
    'n':'index of refraction',
    'n1':'index of refraction1',
    'n2':'index of refraction2',
    'p':'momentum',
    'phi':'work function',
    'theta':'angle',
    'theta1':'angle1',
    'theta2':'angle2',
    'v':'speed'
  }
]

#constant values dict
cons_dict={
    'c':299792458,
    'g':9.81,
    'eo':8.85418782e-12,
    'uo':4e-7 * pi,
    'ke':8.99e9,
    'G':6.6741e-11,
    'hc':6.6261e-34,
    'kb':1.38065e-23,
    'Rc':8.31446
}

#dictionary of symbols for pretty print equation
pps_dict={
  'kappa':'κ',
  'sigma':'σ',
  'theta':'θ',
  'theta1':'θ₁',
  'theta2':'θ₂',
  'flux':'𝚽',
  'sqrt':'√',
  'wave':'λ',
  'phi':'𝜙',
  'tau':'τ',
  '**2':'²',
  'A1':'A₁',
  'A2':'A₂',
  'aa':'α',
  'eo':'ε₀',
  'di':'dᵢ',
  'do':'d₀',
  'h1':'h₁',
  'h2':'h₂',
  'hc':'ℎ',
  'hi':'hᵢ',
  'ho':'h₀',
  'kb':'𝑘B',
  'ke':'kₑ',
  'm1':'m₁',
  'm2':'m₂',
  'n1':'n₁',
  'n2':'n₂',
  'Pi':'Pᵢ',
  'pi':'π',
  'p1':'p₁',
  'p2':'p₂',
  'q1':'q₁',
  'q2':'q₂',
  'Rc':'𝑅',
  'uo':'μ₀',
  'vf':'𝑣𝟋',
  'v1':'𝑣₁',
  'v2':'𝑣₂',
  'wf':'ω𝟋',
  'd':'𝑑',
  'l':'ℓ',
  'p':'ρ',
  'u':'μ',
  'v':'𝑣',
  'w':'ω',
  '*':'⋅'
}

#Equation object
class Eq(object):
  def __init__(self, equation, lst):
    self.var_lst = sorted(lst)
    self.equation = equation
    
  def sol(self,value1,ans_var):
    cprint(ppeq_lst[value1],'yellow')
    ans = solve(self.equation, vl[sect_num][ans_var])
    ans_lst = str(ans).strip("[]")
    ans_lst = ans_lst.split(',')
    #used for chain function
    global temp_ans_lst
    temp_ans_lst = []
    #^used for chain function
    ppans(ans_lst,ans_var)###not implemented print only equations used to get final answer(csolve)
  
  #assists bsolve() to find this equation even with extraneous variables
  def asolve(self,lst,ans_var):
    if set(self.var_lst).issubset(lst) and (ans_var in self.var_lst):
      return True
  
  #odd one out
  def ooo(self,lst):
    temp_lst = list(self.var_lst)
    for var1 in self.var_lst:
      if var1 in lst:
        temp_lst.remove(var1)
    if len(temp_lst) == 1:
      return str(temp_lst).strip("[']")
    else:
      return 0

#section selection help function
def halp():
  print()
  cprint('Here are the available sections of physics: ','magenta')
  for value1, key1 in enumerate(sects):
    if value1 !=  0 and value1 != len(sects):
      print(sects[key1][1] + ': ' + key1)
  print()

#variable selection help function
def vhalp():
  print()
  cprint('Here are the available variables in the ' + sects[sect][1] + ' section:','magenta')
  for key1 in vl3[sect_num]:
    print(key1 + ": " + str(vl3[sect_num][key1]))
  print()
  print()

#prints all equations in sect_num
def print_eq():
  sects[sect][0]('print_eq')

#makes pretty print variables and main() makes ppeq using ppvars, which Eq objects can call
def create_ppeq():
  for key1 in vl[0]:
    vl[0][key1] = Symbol(key1)
  globals().update(vl[0])
  for key1 in vl[sect_num]:
    vl[sect_num][key1] = Symbol(key1)
  globals().update(vl[sect_num])
  sects[sect][0]('create_ppeq')

#custom pretty print replace function
def ppr(string,t1,t2):
  string = string.replace(t1,t2)
  return string

#pretty prints answers with 4 significant figures
def ppans(ans_lst,ans_var):
  for key1 in pps_dict:
    if ans_var == key1:
      ans_var = pps_dict[key1]
  for value1, ans1 in enumerate(ans_lst):
    ans1 = str(sympify(ans1).evalf())
    temp_ans_lst.append(ans1)
    if value1 == len(ans_lst) - 1:
      if '0.00' in ans1:
        print(ans_var + ' = ' + '%.3E' % Decimal(ans1))
      else:
        print(ans_var + ' = ' + str("{0:.4g}".format(float(ans1))))
    else:
      if '0.00' in ans1:
        print(ans_var + ' = ' + '%.3E' % Decimal(ans1), end = ' or ')
      else:
        print(ans_var + ' = ' + str("{0:.4g}".format(float(ans1))), end = ' or ')
  print()  

#creates a set containg subsets of a list that have certain amount of terms
def subify(lst,num):
  return list(combinations(lst,num))

#basic solve function
def bsolve(lst1,lst2):
  score = 0
  
  for value, eq in enumerate(lst1):
    if lst2 == eq.var_lst:
      eq.sol(value,z)
      return 1
    elif eq.asolve(lst2,z):
      eq.sol(value,z)
      return 1 
  return score

#chain solve function
def csolve(lst1,lst2):
  score = 0
  solved = 0
  try:
    lst2.remove(z)
  except ValueError:
    pass
  list4 =[]
  for num in range(1,len(lst2) + 1):
    list4.append(num)
  list4.reverse()
  for num in list4:
    subs = subify(lst2,num)
    for lst3 in subs:
      for value,eq in enumerate(lst1):
        if not set(eq.var_lst).issubset(lst2):
          new_z = eq.ooo(lst3)
          if new_z != 0:
            solved = csolve2(lst1,lst2,value,eq,new_z)
            if solved == 'solved':
              return 'stop'
            else:
              score += solved
      if solved == 'solved':
        return 'stop'
    if solved == 'solved':
      return 'stop'
  return score

#chain solve function part two (once a potential chain equation is found)
def csolve2(lst1,lst2,value,eq,new_z):
  solved = 0
  eq.sol(value,new_z)
  lst2.append(new_z)
  lst2 = sorted(lst2)
  lst1.remove(eq)
  del ppeq_lst[value]
  for ans1 in temp_ans_lst:
    new_z_symbol = Symbol(str(vl[sect_num][new_z]))
    vl[sect_num][new_z] = float(ans1)
    globals().update(vl[sect_num])
    for eq2 in lst1:
      eq2.equation = eq2.equation.subs(new_z_symbol,ans1)
    lst2.append(z)
    lst2 = sorted(lst2)
    solved = sects[sect][0](lst2,lst1)
    if solved == 'solved':
      pass
    elif solved == 'none':
      return 'solved'
  if solved == 'solved':
    return solved
  else:
    return 1

#solver
def solver(lst1,lst2):
  while True:
    score = bsolve(lst1,lst2)
    scorec = 0
    #Chains other equations
    if score != 0:
      return 'solved'
    else:
      scorec = csolve(lst1,lst2)
    
    if scorec == 0:
      cprint('No equation in the database matches your variable inputs.','red',attrs = ['bold'])
      return 'none'
    elif scorec == 'stop':
      return 'none'
    else:
      break

#solve/equation printer/ppeq creator function called at end of each section of physics
def main(lst1, lst2):
  if lst2 == 'print_eq':
    print()
    cprint('Here are the available formulas in the ' + sects[sect][1] + ' section:','magenta')
    for eq in ppeq_lst:
      print(eq)
    print()
  elif lst2 == 'create_ppeq':
    for obj in lst1:
      ppeq = str(obj.equation)
      for key1 in pps_dict:
        ppeq = ppr(ppeq,key1,pps_dict[key1])
      ppeq = '0 = ' + ppeq
      ppeq_lst.append(ppeq)
  #solve equation
  else:
    return solver(lst1, lst2)###





#ALL SECTIONS OF PHYSICS
#--------------------------------
#MECHANICS 
def mech(lst,clist = []):
  eq1 = Eq(-d + (vf + v) / 2 * t, ['d', 't', 'vf', 'v'])
  eq2 = Eq(-vf + v + a * t, ['a', 't', 'vf', 'v'])
  eq3 = Eq(-vf ** 2 + v ** 2 + 2 * a * d, ['a', 'd', 'vf', 'v'])
  eq4 = Eq(-d + v * t + .5 * a * t ** 2, ['a', 'd', 't', 'v'])
  eq5 = Eq(-d + vf * t - .5 * a * t ** 2, ['a', 'd', 't', 'vf'])
  eq6 = Eq(-F + m * a, ['a', 'F', 'm'])
  eq7 = Eq(-a + v ** 2 / r, ['v', 'a', 'r'])
  eq8 = Eq(-p + m * v, ['m','v','p'])
  eq9 = Eq(-p + F * t, ['F','t','p'])
  eq10 = Eq(-K + .5 * m * v, ['m','v','K'])
  eq11 = Eq(-W + F * d, ['W','F','d'])
  eq12 = Eq(-E + F * d, ['E','F','d'])
  eq13 = Eq(-P + W / t, ['P','W','t'])
  eq14 = Eq(-P + E / t, ['P','E','t'])
  eq15 = Eq(-theta + w * t + .5 * aa * t ** 2,['theta','w','t','aa'])
  eq16 = Eq(-wf + w + aa * t,['wf','w','t','aa'])
  eq17 = Eq(-theta + (wf + w) / 2 * t,['theta','wf','w','t'])
  eq18 = Eq(-wf ** 2 + w ** 2 + 2 * aa * theta,['theta','wf','w','aa'])
  eq19 = Eq(-theta + wf * t - .5 * aa * t ** 2,['theta','wf','t','aa'])
  eq20 = Eq(-x + A * cos(w * t),['x','A','w','t'])
  eq21 = Eq(-x + A * cos(2 * pi * f * t),['x','A','f','t'])
  eq22 = Eq(-aa + tau / I,['aa','tau','I'])
  eq23 = Eq(-tau + r * F,['tau','r','F'])
  eq24 = Eq(-tau + r * F * sin(theta),['tau','r','F','theta'])
  eq25 = Eq(-L + I * w,['L','I','w'])
  eq26 = Eq(-L + tau * t,['L','tau','t'])
  eq27 = Eq(-K + .5 * I * w ** 2,['K','I','w'])
  eq28 = Eq(-F + k * x,['F','k','x'])
  eq29 = Eq(-U + .5 * k * x ** 2,['U','k','x'])
  eq30 = Eq(-U + m * g * y,['U','m','g','y'])
  eq31 = Eq(-T + 2 * pi / w,['T','w'])
  eq32 = Eq(-T + 1 / f,['T','f'])
  eq33 = Eq(-w + 2 * pi * f,['w','f'])
  eq34 = Eq(-T + 2 * pi * sqrt(m / k),['T','m','k'])
  eq35 = Eq(-T + 2 * pi * sqrt(l / g),['T','l'])
  eq36 = Eq(-T + 2 * pi * sqrt(l / a),['T','l','a'])
  eq37 = Eq(-F + G * m1 * m2 / r**2,['F','m1','m2','r'])
  eq38 = Eq(-U + G * m1 * m2 / r,['U','m1','m2','r'])
  
  #list of equation objects in mechanics
  mech_eq = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15, eq16, eq17, eq18, eq19, eq20, eq21, eq22, eq23, eq24, eq25, eq26, eq27, eq28, eq29, eq30, eq31, eq32, eq33, eq34, eq35, eq36, eq37, eq38]
  if clist == []:
    main(mech_eq,lst)
  else:
    #csolve is calling this function and has a modified mech_eq list
    return main(clist,lst)

#ELECTRICITY AND MAGNETISM
def em(lst,clist = []):
  eq1 = Eq(-F + ke * q1 * q2 / r ** 2,['F','q1','q2','r'])
  eq2 = Eq(-E + F / q,['E','F','q'])
  eq3 = Eq(-E + ke * q / r ** 2,['E','q','r'])
  eq4 = Eq(-U + q * V,['U','q','V'])
  eq5 = Eq(-V + ke * q / r,['V','q','r'])
  eq6 = Eq(-E + V / d,['E','V','d'])
  eq7 = Eq(-V + Q / C,['V','Q','C'])
  eq8 = Eq(-C + eo * A / d,['C','A','d'])
  eq9 = Eq(-C + kappa * eo * A / d,['C','kappa','A','d'])
  eq10 = Eq(-E + Q / eo / A,['E','Q','A'])
  eq11 = Eq(-U + .5 * C * V ** 2,['U','C','V'])
  eq12 = Eq(-U + .5 * Q * V,['U','Q','V'])
  eq13 = Eq(-U + Q ** 2 / (2 * C),['U','Q','C'])
  eq14 = Eq(-I + Q / t,['I','Q','t'])
  eq15 = Eq(-R + p * l / A,['R','p','l','A'])
  eq16 = Eq(-P + I * V,['P','I','V'])
  eq17 = Eq(-V + I * R,['I','V','R'])
  eq18 = Eq(-B + uo / (2 * pi * r),['B','r'])
  eq19 = Eq(-F + q * v * B,['F','q','v','B'])
  eq20 = Eq(-F + q * v * B * sin(theta),['F','q','v','B','theta'])
  eq21 = Eq(-F + B * I * l,['F','B','I','l'])
  eq22 = Eq(-F + B * I * l * sin(theta),['F','B','I','l','theta'])
  eq23 = Eq(-flux + B * A,['flux','B','A'])
  eq24 = Eq(-flux + B * cos(theta) * A,['flux','B','A','theta'])
  eq25 = Eq(-V - flux / t,['V','flux','t'])
  eq26 = Eq(-V + B * l * v,['V','B','l','v'])
  eq27 = Eq(-U + ke * q1 * q2 / r,['U','q1','q2','r'])
  eq28 = Eq(-v + E / B,['v','E','B'])
  
  #list of equation objects in electricity and magnetism
  em_eq = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15, eq16, eq17, eq18, eq19, eq20, eq21, eq22, eq23, eq24, eq25, eq26, eq27, eq28]
  
  if clist == []:
    main(em_eq,lst)
  else:
    #csolve is calling this function and has a modified em_eq list
    return main(clist,lst)

#FLUID MECHANICS AND THERMAL PHYSICS
def fluid_thermal(lst,clist = []):
  eq1 = Eq(-p + m / V,['p','m','V'])
  eq2 = Eq(-P + F / A,['P','F','A'])
  eq3 = Eq(-P + Pi + p * g * h,['P','Pi','p','h'])
  eq4 = Eq(-P + p * g * h,['P','p','h'])
  eq5 = Eq(-F + p * V * g,['F','p','V',])
  eq6 = Eq(-(A1 * v1) + (A2 * v2),['A1','v1','A2','v2'])
  eq7 = Eq(-(P1 + p * g * h1 + .5 * p * v1 **2) + P2 + p * g * h2 + .5 * p  * v2 ** 2,['P1','P2','p','h1','h2','v1','v2'])
  eq8 = Eq(-(Q / t) + k * A * T / L,['Q','t','k','A','T','L'])
  eq9 = Eq(-(P * V) + n * Rc * T,['P','V','n','T'])
  eq10 = Eq(-(P * V) + N * kb * T,['P','V','N','T'])
  eq11 = Eq(-K + 3 / 2 * kb * T,['K','kb','T'])
  eq12 = Eq(-W - P * V,['W','P','V'])
  eq13 = Eq(-U - Q + W,['U','Q','W'])
  
  ft_eq = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13]
  if clist == []:
    main(ft_eq,lst)
  else:
    #csolve is calling this function and has a modified ft_eq list
    return main(clist,lst)

#MODERN PHYSICS, WAVES, AND OPTICS
def mod(lst,clist = []):
  eq1 = Eq(-wave + v / f,['wave','v','f'])
  eq2 = Eq(-n + c / v,['n','v'])
  eq3 = Eq(-(n1 * sin(theta1)) + n2 * sin(theta2),['n1','n2','theta1','theta2'])
  eq4 = Eq(-1 / f + 1 / di + 1 / do,['f','di','do'])
  eq5 = Eq(-M + hi / ho,['M','hi','ho'])
  eq6 = Eq(-M + di / do,['M','di','do'])
  eq7 = Eq(-L + N * wave,['L','N','wave'])
  eq8 = Eq(-(d * sin(theta)) + N * wave,['d','theta','N','wave'])
  eq9 = Eq(-E + hc * f,['E','f'])
  eq10 = Eq(-K + hc * f - phi,['K','f','phi'])
  eq11 = Eq(-wave + hc / p,['wave','p'])
  eq12 = Eq(-E + m * c ** 2,['E','m'])
  
  mod_eq = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12]
  
  if clist == []:
    main(mod_eq,lst)
  else:
    #csolve is calling this function and has a modified mod_eq list
    return main(clist,lst)





#DICTIONARY OF ALL THE SECTIONS
#--------------------------------
sects={
  '/help':[halp,'HELP'],
  'mech':[mech,'MECHANICS'],
  'e/m':[em,'ELECTRICITY AND MAGNETISM'],
  'fluid':[fluid_thermal,'FLUID MECHANICS AND THERMAL PHYSICS'],
  'mod':[mod,'MODERN PHYSICS, WAVES, AND OPTICS']
}
#--------------------------------





bye = 0
while bye == 0:
  #asking for the section of physics
  def ask_sect():
    global sect
    global sect_num
    global ppeq_lst
    global vl2
    global vl3
    vl3 = []
    while True:
      sect1 = input(colored('Select a section of physics.\n(type /help for help, /alleq to display all equations, /allvar to display all variables)\n','cyan'))
      sect1 = sect1.lower()
      for sect2 in sects:
        sect_num = list(sects.keys()).index(sect2)
        vl2 = dict(vl[sect_num])
        vl3.append(dict(vl2))
      if sect1 == '/help':
        sects[sect1][0]()
      elif sect1 == '/alleq':
        for sect2 in sects:
          if sect2 != '/help':
            sect = sect2
            sect_num = list(sects.keys()).index(sect2)
            ppeq_lst = []
            create_ppeq()
            print_eq()
      elif sect1 == '/allvar':
        for sect2 in sects:
          if sect2 != '/help':
            sect = sect2
            sect_num = list(sects.keys()).index(sect2)
            vhalp()
      elif sect1 in sects:
        return sect1
      else:
        print('That is not an available section.')
  sect = ask_sect()
  
  #sect_num links vl and sects both have their index 0 as something else so it works fine
  sect_num = list(sects.keys()).index(sect)
  
  #creating a back up vl to print if halp() is called
  vl2 = dict(vl[sect_num])
  
  #creating list that will contain pretty print equations for vl[sect_num]
  ppeq_lst = []
  
  #creating pretty print equations
  create_ppeq()
  
  #asking for known variables
  def ask_vars(sect_num1):
    while True:
      y1 = input(colored('what variables are known?\n(Enter with commas in between)(type /var to list variables, /eq to display equations)\n','cyan'))
      y1 = sorted([var.strip() for var in y1.split(',')])
      score1 = 0
      for var1 in y1:
        if var1 in vl[sect_num1]:
          score1 +=1
      if y1 == ['/var']:
        vhalp()
      elif y1 == ['/eq']:
        print_eq()
      elif score1 == len(y1):
        return y1
      else:
        cprint('One or more of those is not an available variable.','red',attrs = ['bold'])
  y_lst = ask_vars(sect_num)
  
  #asking for the unknown variable
  def ask_var(sect_num1):
    while True:
      z1 = input(colored('what variable are you solving for?\n(type /var to list variables, /eq to display equations)\n','cyan'))
      z1 = z1.strip()
      
      if z1 == '/var':
        vhalp()
      elif z1 == '/eq':
        print_eq()
      elif z1 in y_lst:
        print('You inputed that as a known variable.')
      elif z1 in vl[sect_num1]:
        return z1
      else:
        cprint("That is not an available variable.",'red',attrs = ['bold'])
  z = ask_var(sect_num)
  
  #creating a list containing known and unknown variables
  yz = list(y_lst)
  yz.append(z)
  yz = sorted(yz)
  
  #creating a 2nd yz to test if csolve was used to update yz list
  yz2 = list(yz)
  
  #updating all values (not keys) in vl[sect] into sympy symbols to prevent error
  for key in vl[sect_num]:
    vl[sect_num][key] = Symbol(key)
  
  #updating all constant values in vl[0] into their respective values
  for key in vl[0]:
    vl[0][key] = cons_dict[key]
  
  #assigning floats to variables
  cprint('Please enter the values of your variables:','cyan')
  def assign(ele1):
    while True:
      try:
        return float(input(str(ele1) + ' ='))
      except ValueError:
        print("Please enter a number.")
  for ele in y_lst:
    ppele = str(ele)
    for key in pps_dict:
      if ele == key:
        ppele = pps_dict[key]
    vl[sect_num][ele] = assign(ppele)
  
  #store the variable we are solving for in vl[z] and set as x in sympy
  vl[sect_num][z] = Symbol(z)
  
  #updating the keys in dicitonary vl into variables (pretty cool)
  error = 'who is the physics god?'
  globals().update(vl[0])
  globals().update(vl[sect_num])
  
  #solve
  sects[sect][0](yz)
  
  again = input('solve another equation? (yes/no)\n')
  if again == ('no' or 'No'):
    bye += 1
  elif again == ('yes' or 'Yes'):
    pass
  elif again == error:
    cprint('ALBERT','red',attrs=['bold','blink','underline'])
  else:
    cprint("c'mon, u can't even type yes or no? Goodbye.",'red',attrs = ['bold'])
    bye += 1
