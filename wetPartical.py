import numpy as np 

c_ss = {# sea salt
    "c1":0.7674, 
    "c2":3.079,
    "c3":2.573e-11,
    "c4":-1.424
}
c_ur = {# urban
    "c1":0.3926, 
    "c2":3.101,
    "c3":4.190e-11,
    "c4":-1.404
}
c_ru = {# rural
    "c1":0.2789, 
    "c2":3.115,
    "c3":5.415e-11,
    "c4":-1.399
}
c_ns = {# (NH4)2SO4
    "c1":0.4809, 
    "c2":3.082,
    "c3":3.110e-11,
    "c4":-1.428
}
pm_type_dict = {
    'ss':c_ss,
    'ur':c_ur,
    'ru':c_ru,
    'ns':c_ns
}

def calcu_Dw(rh,Dd,pm_type):
    '''计算湿颗粒直径 (Gerber, 1985)
    rh:相对湿度
    Dd:干颗粒直径
    pm_type:颗粒类型，分为：
        ss：sea salt
        ur：urban
        ru：rural
        ns：(NH4)2SO4
    '''    
    var = pm_type_dict[pm_type]
    c1, c2, c3, c4 = var['c1'], var['c2'], var['c3'], var['c4']
    Dw = 2 * (c1*(Dd/2)**c2/(c3*(Dd/2)**c4-np.log10(rh)) + (Dd/2)**3) **(1/3)
    return Dw
print(calcu_Dw(rh=0.68, Dd=2.5, pm_type='ss'))