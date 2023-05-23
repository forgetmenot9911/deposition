from deposition import calc_Dw_Gerber1985,calc_Dw_Quinn1998,\
                        calc_vd,calc_ustar,cc,\
                        c_ss, c_ns, c_ru, c_ur, pm_type_dict,\
                        wspd, z, Ta, Tw, rh, u_starIni, rho, \
                            Dd, Kc, k, g, err, nu, phi 
import pandas as pd 

def main():
    df = pd.read_excel('Particle_Sizes_of_Aerosol_Iron.xlsx',
                       header=2)

if __name__ == '__main__':
    main()