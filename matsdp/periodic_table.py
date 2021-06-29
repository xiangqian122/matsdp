# -*- coding: utf-8 -*-
def periodic_tab():
  periodic_table_dict = {}
  periodic_table_dict['atomic_number']={'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}
  periodic_table_dict['symbol']={'H': 'H', 'He': 'He', 'Li': 'Li', 'Be': 'Be', 'B': 'B', 'C': 'C', 'N': 'N', 'O': 'O', 'F': 'F', 'Ne': 'Ne', 'Na': 'Na', 'Mg': 'Mg', 'Al': 'Al', 'Si': 'Si', 'P': 'P', 'S': 'S', 'Cl': 'Cl', 'Ar': 'Ar', 'K': 'K', 'Ca': 'Ca', 'Sc': 'Sc', 'Ti': 'Ti', 'V': 'V', 'Cr': 'Cr', 'Mn': 'Mn', 'Fe': 'Fe', 'Co': 'Co', 'Ni': 'Ni', 'Cu': 'Cu', 'Zn': 'Zn', 'Ga': 'Ga', 'Ge': 'Ge', 'As': 'As', 'Se': 'Se', 'Br': 'Br', 'Kr': 'Kr', 'Rb': 'Rb', 'Sr': 'Sr', 'Y': 'Y', 'Zr': 'Zr', 'Nb': 'Nb', 'Mo': 'Mo', 'Tc': 'Tc', 'Ru': 'Ru', 'Rh': 'Rh', 'Pd': 'Pd', 'Ag': 'Ag', 'Cd': 'Cd', 'In': 'In', 'Sn': 'Sn', 'Sb': 'Sb', 'Te': 'Te', 'I': 'I', 'Xe': 'Xe', 'Cs': 'Cs', 'Ba': 'Ba', 'La': 'La', 'Ce': 'Ce', 'Pr': 'Pr', 'Nd': 'Nd', 'Pm': 'Pm', 'Sm': 'Sm', 'Eu': 'Eu', 'Gd': 'Gd', 'Tb': 'Tb', 'Dy': 'Dy', 'Ho': 'Ho', 'Er': 'Er', 'Tm': 'Tm', 'Yb': 'Yb', 'Lu': 'Lu', 'Hf': 'Hf', 'Ta': 'Ta', 'W': 'W', 'Re': 'Re', 'Os': 'Os', 'Ir': 'Ir', 'Pt': 'Pt', 'Au': 'Au', 'Hg': 'Hg', 'Tl': 'Tl', 'Pb': 'Pb', 'Bi': 'Bi', 'Po': 'Po', 'At': 'At', 'Rn': 'Rn', 'Fr': 'Fr', 'Ra': 'Ra', 'Ac': 'Ac', 'Th': 'Th', 'Pa': 'Pa', 'U': 'U', 'Np': 'Np', 'Pu': 'Pu', 'Am': 'Am', 'Cm': 'Cm', 'Bk': 'Bk', 'Cf': 'Cf', 'Es': 'Es', 'Fm': 'Fm', 'Md': 'Md', 'No': 'No', 'Lr': 'Lr', 'Rf': 'Rf', 'Db': 'Db', 'Sg': 'Sg', 'Bh': 'Bh', 'Hs': 'Hs', 'Mt': 'Mt', 'Ds': 'Ds', 'Rg': 'Rg', 'Cn': 'Cn', 'Nh': 'Nh', 'Fl': 'Fl', 'Mc': 'Mc', 'Lv': 'Lv', 'Ts': 'Ts', 'Og': 'Og'}
  periodic_table_dict['name']={'H': 'Hydrogen', 'He': 'Helium', 'Li': 'Lithium', 'Be': 'Beryllium', 'B': 'Boron', 'C': 'Carbon', 'N': 'Nitrogen', 'O': 'Oxygen', 'F': 'Fluorine', 'Ne': 'Neon', 'Na': 'Sodium', 'Mg': 'Magnesium', 'Al': 'Aluminum', 'Si': 'Silicon', 'P': 'Phosphorus', 'S': 'Sulfur', 'Cl': 'Chlorine', 'Ar': 'Argon', 'K': 'Potassium', 'Ca': 'Calcium', 'Sc': 'Scandium', 'Ti': 'Titanium', 'V': 'Vanadium', 'Cr': 'Chromium', 'Mn': 'Manganese', 'Fe': 'Iron', 'Co': 'Cobalt', 'Ni': 'Nickel', 'Cu': 'Copper', 'Zn': 'Zinc', 'Ga': 'Gallium', 'Ge': 'Germanium', 'As': 'Arsenic', 'Se': 'Selenium', 'Br': 'Bromine', 'Kr': 'Krypton', 'Rb': 'Rubidium', 'Sr': 'Strontium', 'Y': 'Yttrium', 'Zr': 'Zirconium', 'Nb': 'Niobium', 'Mo': 'Molybdenum', 'Tc': 'Technetium', 'Ru': 'Ruthenium', 'Rh': 'Rhodium', 'Pd': 'Palladium', 'Ag': 'Silver', 'Cd': 'Cadmium', 'In': 'Indium', 'Sn': 'Tin', 'Sb': 'Antimony', 'Te': 'Tellurium', 'I': 'Iodine', 'Xe': 'Xenon', 'Cs': 'Cesium', 'Ba': 'Barium', 'La': 'Lanthanum', 'Ce': 'Cerium', 'Pr': 'Praseodymium', 'Nd': 'Neodymium', 'Pm': 'Promethium', 'Sm': 'Samarium', 'Eu': 'Europium', 'Gd': 'Gadolinium', 'Tb': 'Terbium', 'Dy': 'Dysprosium', 'Ho': 'Holmium', 'Er': 'Erbium', 'Tm': 'Thulium', 'Yb': 'Ytterbium', 'Lu': 'Lutetium', 'Hf': 'Hafnium', 'Ta': 'Tantalum', 'W': 'Tungsten', 'Re': 'Rhenium', 'Os': 'Osmium', 'Ir': 'Iridium', 'Pt': 'Platinum', 'Au': 'Gold', 'Hg': 'Mercury', 'Tl': 'Thallium', 'Pb': 'Lead', 'Bi': 'Bismuth', 'Po': 'Polonium', 'At': 'Astatine', 'Rn': 'Radon', 'Fr': 'Francium', 'Ra': 'Radium', 'Ac': 'Actinium', 'Th': 'Thorium', 'Pa': 'Protactinium', 'U': 'Uranium', 'Np': 'Neptunium', 'Pu': 'Plutonium', 'Am': 'Americium', 'Cm': 'Curium', 'Bk': 'Berkelium', 'Cf': 'Californium', 'Es': 'Einsteinium', 'Fm': 'Fermium', 'Md': 'Mendelevium', 'No': 'Nobelium', 'Lr': 'Lawrencium', 'Rf': 'Rutherfordium', 'Db': 'Dubnium', 'Sg': 'Seaborgium', 'Bh': 'Bohrium', 'Hs': 'Hassium', 'Mt': 'Meitnerium', 'Ds': 'Darmstadtium', 'Rg': 'Roentgenium', 'Cn': 'Copernicium', 'Nh': 'Nihonium', 'Fl': 'Flerovium', 'Mc': 'Moscovium', 'Lv': 'Livermorium', 'Ts': 'Tennessine', 'Og': 'Oganesson'}
  periodic_table_dict['elmt_color']={'H': 'gray', 'He': 'gray', 'Li': 'gray', 'Be': 'gray', 'B': 'gray', 'C': 'gray', 'N': 'gray', 'O': 'gray', 'F': 'gray', 'Ne': 'gray', 'Na': 'gray', 'Mg': 'gray', 'Al': 'black', 'Si': 'gray', 'P': 'gray', 'S': 'gray', 'Cl': 'gray', 'Ar': 'gray', 'K': 'gray', 'Ca': 'gray', 'Sc': 'lemonchiffon', 'Ti': 'rosybrown', 'V': 'wheat', 'Cr': 'lightgreen', 'Mn': 'chocolate', 'Fe': 'lightgray', 'Co': 'purple', 'Ni': 'gray', 'Cu': 'peru', 'Zn': 'turquoise', 'Ga': 'gray', 'Ge': 'gray', 'As': 'gray', 'Se': 'gray', 'Br': 'gray', 'Kr': 'gray', 'Rb': 'gray', 'Sr': 'gray', 'Y': 'springgreen', 'Zr': 'sienna', 'Nb': 'firebrick', 'Mo': 'lightskyblue', 'Tc': 'lawngreen', 'Ru': 'cyan', 'Rh': 'darkslateblue', 'Pd': 'orchid', 'Ag': 'silver', 'Cd': 'olive', 'In': 'gray', 'Sn': 'gray', 'Sb': 'gray', 'Te': 'gray', 'I': 'gray', 'Xe': 'gray', 'Cs': 'gray', 'Ba': 'gray', 'La': 'gray', 'Ce': 'gray', 'Pr': 'gray', 'Nd': 'gray', 'Pm': 'gray', 'Sm': 'gray', 'Eu': 'gray', 'Gd': 'gray', 'Tb': 'gray', 'Dy': 'gray', 'Ho': 'gray', 'Er': 'gray', 'Tm': 'gray', 'Yb': 'gray', 'Lu': 'gray', 'Hf': 'yellow', 'Ta': 'magenta', 'W': 'blue', 'Re': 'red', 'Os': 'brown', 'Ir': 'plum', 'Pt': 'aqua', 'Au': 'gold', 'Hg': 'mediumblue', 'Tl': 'gray', 'Pb': 'gray', 'Bi': 'gray', 'Po': 'gray', 'At': 'gray', 'Rn': 'gray', 'Fr': 'gray', 'Ra': 'gray', 'Ac': 'gray', 'Th': 'gray', 'Pa': 'gray', 'U': 'gray', 'Np': 'gray', 'Pu': 'gray', 'Am': 'gray', 'Cm': 'gray', 'Bk': 'gray', 'Cf': 'gray', 'Es': 'gray', 'Fm': 'gray', 'Md': 'gray', 'No': 'gray', 'Lr': 'gray', 'Rf': 'gray', 'Db': 'gray', 'Sg': 'gray', 'Bh': 'gray', 'Hs': 'gray', 'Mt': 'gray', 'Ds': 'gray', 'Rg': 'gray', 'Cn': 'gray', 'Nh': 'gray', 'Fl': 'gray', 'Mc': 'gray', 'Lv': 'gray', 'Ts': 'gray', 'Og': 'gray'}
  periodic_table_dict['atomic_radius']={'H': 25.0, 'He': 31.0, 'Li': 145.0, 'Be': 105.0, 'B': 85.0, 'C': 70.0, 'N': 65.0, 'O': 60.0, 'F': 50.0, 'Ne': 51.0, 'Na': 180.0, 'Mg': 150.0, 'Al': 125.0, 'Si': 110.0, 'P': 100.0, 'S': 100.0, 'Cl': 100.0, 'Ar': 71.0, 'K': 220.0, 'Ca': 180.0, 'Sc': 160.0, 'Ti': 140.0, 'V': 135.0, 'Cr': 140.0, 'Mn': 140.0, 'Fe': 140.0, 'Co': 135.0, 'Ni': 135.0, 'Cu': 135.0, 'Zn': 135.0, 'Ga': 130.0, 'Ge': 125.0, 'As': 115.0, 'Se': 115.0, 'Br': 115.0, 'Kr': 103.0, 'Rb': 235.0, 'Sr': 200.0, 'Y': 180.0, 'Zr': 155.0, 'Nb': 145.0, 'Mo': 145.0, 'Tc': 135.0, 'Ru': 130.0, 'Rh': 135.0, 'Pd': 140.0, 'Ag': 160.0, 'Cd': 155.0, 'In': 155.0, 'Sn': 145.0, 'Sb': 145.0, 'Te': 140.0, 'I': 140.0, 'Xe': 124.0, 'Cs': 260.0, 'Ba': 215.0, 'La': 195.0, 'Ce': 185.0, 'Pr': 185.0, 'Nd': 185.0, 'Pm': 185.0, 'Sm': 185.0, 'Eu': 185.0, 'Gd': 180.0, 'Tb': 175.0, 'Dy': 175.0, 'Ho': 175.0, 'Er': 175.0, 'Tm': 175.0, 'Yb': 175.0, 'Lu': 175.0, 'Hf': 155.0, 'Ta': 145.0, 'W': 135.0, 'Re': 135.0, 'Os': 130.0, 'Ir': 135.0, 'Pt': 135.0, 'Au': 135.0, 'Hg': 150.0, 'Tl': 190.0, 'Pb': 180.0, 'Bi': 160.0, 'Po': 190.0, 'At': None, 'Rn': 134.0, 'Fr': 260.0, 'Ra': 215.0, 'Ac': 195.0, 'Th': 180.0, 'Pa': 180.0, 'U': 175.0, 'Np': 175.0, 'Pu': 175.0, 'Am': 175.0, 'Cm': 135.0, 'Bk': None, 'Cf': None, 'Es': None, 'Fm': None, 'Md': None, 'No': None, 'Lr': None, 'Rf': None, 'Db': None, 'Sg': None, 'Bh': None, 'Hs': None, 'Mt': None, 'Ds': None, 'Rg': None, 'Cn': None, 'Nh': None, 'Fl': None, 'Mc': None, 'Lv': None, 'Ts': None, 'Og': None}
  periodic_table_dict['metallic_radius']={'H': None, 'He': None, 'Li': 152.0, 'Be': 112.0, 'B': None, 'C': 197.0, 'N': None, 'O': None, 'F': None, 'Ne': None, 'Na': 186.0, 'Mg': 160.0, 'Al': 143.0, 'Si': None, 'P': None, 'S': None, 'Cl': None, 'Ar': None, 'K': 227.0, 'Ca': 197.0, 'Sc': 162.0, 'Ti': 147.0, 'V': 134.0, 'Cr': 128.0, 'Mn': 127.0, 'Fe': 126.0, 'Co': 125.0, 'Ni': 124.0, 'Cu': 128.0, 'Zn': 134.0, 'Ga': 135.0, 'Ge': None, 'As': None, 'Se': None, 'Br': None, 'Kr': None, 'Rb': 248.0, 'Sr': 215.0, 'Y': 180.0, 'Zr': 160.0, 'Nb': 146.0, 'Mo': 139.0, 'Tc': 136.0, 'Ru': 134.0, 'Rh': 134.0, 'Pd': 137.0, 'Ag': 144.0, 'Cd': 151.0, 'In': 167.0, 'Sn': None, 'Sb': None, 'Te': None, 'I': None, 'Xe': None, 'Cs': 265.0, 'Ba': 222.0, 'La': 187.0, 'Ce': 181.8, 'Pr': 182.4, 'Nd': 181.4, 'Pm': 183.4, 'Sm': 180.4, 'Eu': 208.4, 'Gd': 180.4, 'Tb': 177.3, 'Dy': 178.1, 'Ho': 176.2, 'Er': 176.1, 'Tm': 175.9, 'Yb': 193.3, 'Lu': 173.8, 'Hf': 159.0, 'Ta': 146.0, 'W': 139.0, 'Re': 137.0, 'Os': 135.0, 'Ir': 135.0, 'Pt': 138.5, 'Au': 144.0, 'Hg': 151.0, 'Tl': 170.0, 'Pb': None, 'Bi': None, 'Po': None, 'At': None, 'Rn': None, 'Fr': None, 'Ra': None, 'Ac': None, 'Th': 179.0, 'Pa': 263.0, 'U': 165.0, 'Np': 155.0, 'Pu': 159.0, 'Am': 173.0, 'Cm': 174.0, 'Bk': 170.0, 'Cf': 186.0, 'Es': 186.0, 'Fm': None, 'Md': None, 'No': None, 'Lr': None, 'Rf': None, 'Db': None, 'Sg': None, 'Bh': None, 'Hs': 118.0, 'Mt': None, 'Ds': None, 'Rg': None, 'Cn': None, 'Nh': None, 'Fl': None, 'Mc': None, 'Lv': None, 'Ts': None, 'Og': None}
  periodic_table_dict['csd_covalent_radius']={'H': 31.0, 'He': 28.0, 'Li': 128.0, 'Be': 96.0, 'B': 84.0, 'C': 76.0, 'N': 71.0, 'O': 66.0, 'F': 57.0, 'Ne': 58.0, 'Na': 66.0, 'Mg': 141.0, 'Al': 121.0, 'Si': 111.0, 'P': 107.0, 'S': 105.0, 'Cl': 102.0, 'Ar': 106.0, 'K': 203.0, 'Ca': 176.0, 'Sc': 170.0, 'Ti': 160.0, 'V': 153.0, 'Cr': 139.0, 'Mn': 139.0, 'Fe': 132.0, 'Co': 126.0, 'Ni': 124.0, 'Cu': 132.0, 'Zn': 122.0, 'Ga': 122.0, 'Ge': 120.0, 'As': 119.0, 'Se': 120.0, 'Br': 120.0, 'Kr': 116.0, 'Rb': 220.0, 'Sr': 195.0, 'Y': 190.0, 'Zr': 175.0, 'Nb': 164.0, 'Mo': 154.0, 'Tc': 147.0, 'Ru': 146.0, 'Rh': 142.0, 'Pd': 139.0, 'Ag': 145.0, 'Cd': 144.0, 'In': 142.0, 'Sn': 139.0, 'Sb': 139.0, 'Te': 138.0, 'I': 139.0, 'Xe': 140.0, 'Cs': 244.0, 'Ba': 215.0, 'La': 207.0, 'Ce': 204.0, 'Pr': 203.0, 'Nd': 201.0, 'Pm': 199.0, 'Sm': 198.0, 'Eu': 198.0, 'Gd': 196.0, 'Tb': 194.0, 'Dy': 192.0, 'Ho': 192.0, 'Er': 189.0, 'Tm': 190.0, 'Yb': 187.0, 'Lu': 187.0, 'Hf': 175.0, 'Ta': 170.0, 'W': 162.0, 'Re': 151.0, 'Os': 144.0, 'Ir': 141.0, 'Pt': 136.0, 'Au': 136.0, 'Hg': 132.0, 'Tl': 145.0, 'Pb': 146.0, 'Bi': 148.0, 'Po': 140.0, 'At': 150.0, 'Rn': 150.0, 'Fr': 260.0, 'Ra': 221.0, 'Ac': 215.0, 'Th': 206.0, 'Pa': 200.0, 'U': 196.0, 'Np': 190.0, 'Pu': 187.0, 'Am': 180.0, 'Cm': 169.0, 'Bk': None, 'Cf': None, 'Es': None, 'Fm': None, 'Md': None, 'No': None, 'Lr': None, 'Rf': None, 'Db': None, 'Sg': None, 'Bh': None, 'Hs': None, 'Mt': None, 'Ds': None, 'Rg': None, 'Cn': None, 'Nh': None, 'Fl': None, 'Mc': None, 'Lv': None, 'Ts': None, 'Og': None}
  periodic_table_dict['vasppot_paw']={'H': '', 'He': '', 'Li': '_sv', 'Be': '', 'B': '', 'C': '', 'N': '', 'O': '', 'F': '', 'Ne': '', 'Na': '_pv', 'Mg': '', 'Al': '', 'Si': '', 'P': '', 'S': '', 'Cl': '', 'Ar': '', 'K': '_sv', 'Ca': '_sv', 'Sc': '_sv', 'Ti': '_sv', 'V': '', 'Cr': '_pv', 'Mn': '_pv', 'Fe': '', 'Co': '', 'Ni': '', 'Cu': '', 'Zn': '', 'Ga': '_d', 'Ge': '', 'As': '_d', 'Se': '', 'Br': '', 'Kr': '', 'Rb': '_sv', 'Sr': '_sv', 'Y': '_sv', 'Zr': '_sv', 'Nb': '_sv', 'Mo': '_sv', 'Tc': '_sv', 'Ru': '_pv', 'Rh': '_pv', 'Pd': '', 'Ag': '', 'Cd': '', 'In': '_d', 'Sn': '_d', 'Sb': '', 'Te': '', 'I': '', 'Xe': '', 'Cs': '_sv', 'Ba': '_sv', 'La': '', 'Ce': '', 'Pr': '', 'Nd': '', 'Pm': '', 'Sm': '', 'Eu': '', 'Gd': '', 'Tb': '', 'Dy': '', 'Ho': '', 'Er': '', 'Tm': '', 'Yb': '', 'Lu': '', 'Hf': '_pv', 'Ta': '_pv', 'W': '_pv', 'Re': '', 'Os': '', 'Ir': '', 'Pt': '', 'Au': '', 'Hg': '', 'Tl': '_d', 'Pb': '_d', 'Bi': '_d', 'Po': '_d', 'At': '_d', 'Rn': '', 'Fr': '_sv', 'Ra': '_sv', 'Ac': '', 'Th': '', 'Pa': '', 'U': '', 'Np': '', 'Pu': '', 'Am': '', 'Cm': '', 'Bk': '', 'Cf': '', 'Es': '', 'Fm': '', 'Md': '', 'No': '', 'Lr': '', 'Rf': '', 'Db': '', 'Sg': '', 'Bh': '', 'Hs': '', 'Mt': '', 'Ds': '', 'Rg': '', 'Cn': '', 'Nh': '', 'Fl': '', 'Mc': '', 'Lv': '', 'Ts': '', 'Og': ''}
  periodic_table_dict['vasppot_gw']={'H': '_GW', 'He': '_GW', 'Li': '_sv_GW', 'Be': '_sv_GW', 'B': '_GW', 'C': '_GW', 'N': '_GW', 'O': '_GW', 'F': '_GW', 'Ne': '_GW', 'Na': '_sv_GW', 'Mg': '_sv_GW', 'Al': '_GW', 'Si': '_GW', 'P': '_GW', 'S': '_GW', 'Cl': '_GW', 'Ar': '_GW', 'K': '_sv_GW', 'Ca': '_sv_GW', 'Sc': '_sv_GW', 'Ti': '_sv_GW', 'V': '_sv_GW', 'Cr': '_sv_GW', 'Mn': '_sv_GW', 'Fe': '_sv_GW', 'Co': '_sv_GW', 'Ni': '_sv_GW', 'Cu': '_sv_GW', 'Zn': '_sv_GW', 'Ga': '_d_GW', 'Ge': '_d_GW', 'As': '_GW', 'Se': '_GW', 'Br': '_GW', 'Kr': '_GW', 'Rb': '_sv_GW', 'Sr': '_sv_GW', 'Y': '_sv_GW', 'Zr': '_sv_GW', 'Nb': '_sv_GW', 'Mo': '_sv_GW', 'Tc': '_sv_GW', 'Ru': '_sv_GW', 'Rh': '_sv_GW', 'Pd': '_sv_GW', 'Ag': '_sv_GW', 'Cd': '_sv_GW', 'In': '_d_GW', 'Sn': '_d_GW', 'Sb': '_d_GW', 'Te': '_GW', 'I': '_GW', 'Xe': '_GW', 'Cs': '_sv_GW', 'Ba': '_sv_GW', 'La': '_GW', 'Ce': '_GW', 'Pr': '_GW', 'Nd': '_GW', 'Pm': '_GW', 'Sm': '_GW', 'Eu': '_GW', 'Gd': '_GW', 'Tb': '_GW', 'Dy': '_GW', 'Ho': '_GW', 'Er': '_GW', 'Tm': '_GW', 'Yb': '_GW', 'Lu': '_GW', 'Hf': '_sv_GW', 'Ta': '_sv_GW', 'W': '_sv_GW', 'Re': '_sv_GW', 'Os': '_sv_GW', 'Ir': '_sv_GW', 'Pt': '_sv_GW', 'Au': '_sv_GW', 'Hg': '_sv_GW', 'Tl': '_d_GW', 'Pb': '_d_GW', 'Bi': '_d_GW', 'Po': '_d_GW', 'At': '_d_GW', 'Rn': '_d_GW', 'Fr': '_GW', 'Ra': '_GW', 'Ac': '_GW', 'Th': '_GW', 'Pa': '_GW', 'U': '_GW', 'Np': '_GW', 'Pu': '_GW', 'Am': '_GW', 'Cm': '_GW', 'Bk': '_GW', 'Cf': '_GW', 'Es': '_GW', 'Fm': '_GW', 'Md': '_GW', 'No': '_GW', 'Lr': '_GW', 'Rf': '_GW', 'Db': '_GW', 'Sg': '_GW', 'Bh': '_GW', 'Hs': '_GW', 'Mt': '_GW', 'Ds': '_GW', 'Rg': '_GW', 'Cn': '_GW', 'Nh': '_GW', 'Fl': '_GW', 'Mc': '_GW', 'Lv': '_GW', 'Ts': '_GW', 'Og': '_GW'}
  periodic_table_dict['dos_mode']={'H': ['s'], 'He': ['s'], 'Li': ['s'], 'Be': ['s'], 'B': ['s', 'p'], 'C': ['s', 'p'], 'N': ['s', 'p'], 'O': ['s', 'p'], 'F': ['s', 'p'], 'Ne': ['s', 'p'], 'Na': ['s', 'p'], 'Mg': ['s', 'p'], 'Al': ['s', 'p'], 'Si': ['s', 'p'], 'P': ['s', 'p'], 'S': ['s', 'p'], 'Cl': ['s', 'p'], 'Ar': ['s', 'p'], 'K': ['s', 'p'], 'Ca': ['s', 'p'], 'Sc': ['s', 'p', 'd'], 'Ti': ['s', 'p', 'd'], 'V': ['s', 'p', 'd'], 'Cr': ['s', 'p', 'd'], 'Mn': ['s', 'p', 'd'], 'Fe': ['s', 'p', 'd'], 'Co': ['s', 'p', 'd'], 'Ni': ['s', 'p', 'd'], 'Cu': ['s', 'p', 'd'], 'Zn': ['s', 'p', 'd'], 'Ga': ['s', 'p', 'd'], 'Ge': ['s', 'p', 'd'], 'As': ['s', 'p', 'd'], 'Se': ['s', 'p', 'd'], 'Br': ['s', 'p', 'd'], 'Kr': ['s', 'p', 'd'], 'Rb': ['s', 'p', 'd'], 'Sr': ['s', 'p', 'd'], 'Y': ['s', 'p', 'd'], 'Zr': ['s', 'p', 'd'], 'Nb': ['s', 'p', 'd'], 'Mo': ['s', 'p', 'd'], 'Tc': ['s', 'p', 'd'], 'Ru': ['s', 'p', 'd'], 'Rh': ['s', 'p', 'd'], 'Pd': ['s', 'p', 'd'], 'Ag': ['s', 'p', 'd'], 'Cd': ['s', 'p', 'd'], 'In': ['s', 'p', 'd'], 'Sn': ['s', 'p', 'd'], 'Sb': ['s', 'p', 'd'], 'Te': ['s', 'p', 'd'], 'I': ['s', 'p', 'd'], 'Xe': ['s', 'p', 'd'], 'Cs': ['s', 'p', 'd'], 'Ba': ['s', 'p', 'd'], 'La': ['s', 'p', 'd'], 'Ce': ['s', 'p', 'd'], 'Pr': ['s', 'p', 'd'], 'Nd': ['s', 'p', 'd'], 'Pm': ['s', 'p', 'd'], 'Sm': ['s', 'p', 'd'], 'Eu': ['s', 'p', 'd'], 'Gd': ['s', 'p', 'd'], 'Tb': ['s', 'p', 'd'], 'Dy': ['s', 'p', 'd'], 'Ho': ['s', 'p', 'd'], 'Er': ['s', 'p', 'd'], 'Tm': ['s', 'p', 'd'], 'Yb': ['s', 'p', 'd'], 'Lu': ['s', 'p', 'd'], 'Hf': ['s', 'p', 'd'], 'Ta': ['s', 'p', 'd'], 'W': ['s', 'p', 'd'], 'Re': ['s', 'p', 'd'], 'Os': ['s', 'p', 'd'], 'Ir': ['s', 'p', 'd'], 'Pt': ['s', 'p', 'd'], 'Au': ['s', 'p', 'd'], 'Hg': ['s', 'p', 'd'], 'Tl': ['s', 'p', 'd'], 'Pb': ['s', 'p', 'd'], 'Bi': ['s', 'p', 'd'], 'Po': ['s', 'p', 'd'], 'At': ['s', 'p', 'd'], 'Rn': ['s', 'p', 'd'], 'Fr': ['s', 'p', 'd'], 'Ra': ['s', 'p', 'd'], 'Ac': ['s', 'p', 'd'], 'Th': ['s', 'p', 'd'], 'Pa': ['s', 'p', 'd'], 'U': ['s', 'p', 'd'], 'Np': ['s', 'p', 'd'], 'Pu': ['s', 'p', 'd'], 'Am': ['s', 'p', 'd'], 'Cm': ['s', 'p', 'd'], 'Bk': ['s', 'p', 'd'], 'Cf': ['s', 'p', 'd'], 'Es': ['s', 'p', 'd'], 'Fm': ['s', 'p', 'd'], 'Md': ['s', 'p', 'd'], 'No': ['s', 'p', 'd'], 'Lr': ['s', 'p', 'd'], 'Rf': ['s', 'p', 'd'], 'Db': ['s', 'p', 'd'], 'Sg': ['s', 'p', 'd'], 'Bh': ['s', 'p', 'd'], 'Hs': ['s', 'p', 'd'], 'Mt': ['s', 'p', 'd'], 'Ds': ['s', 'p', 'd'], 'Rg': ['s', 'p', 'd'], 'Cn': ['s', 'p', 'd'], 'Nh': ['s', 'p', 'd'], 'Fl': ['s', 'p', 'd'], 'Mc': ['s', 'p', 'd'], 'Lv': ['s', 'p', 'd'], 'Ts': ['s', 'p', 'd'], 'Og': ['s', 'p', 'd']}
  periodic_table_dict['valence_elec']={'H': '1s1', 'He': '1s2', 'Li': '1s1-2s2-2p0', 'Be': '2s2-2p0', 'B': '2s2-2p1', 'C': '2s2-2p2', 'N': '2s2-2p3', 'O': '2s2-2p4', 'F': '2s2-2p5', 'Ne': '2s2-2p6', 'Na': '2p6-3s1', 'Mg': '3s2-3p0', 'Al': '3s2-3p1', 'Si': '3s2-3p2', 'P': '3s2-3p3', 'S': '3s2-3p4', 'Cl': '3s2-3p5', 'Ar': '3s2-3p6', 'K': '3s2-3p6-4s1', 'Ca': '3s2-3p6-4s2', 'Sc': '3s2-3p6-3d1-4s2', 'Ti': '3s2-3p6-4s2-3d2', 'V': '3p6-4d3-4s2', 'Cr': '3p6-3d5-4s1', 'Mn': '3p6-4s2-3d5', 'Fe': '3d6-4s2', 'Co': '3d7-4s2', 'Ni': '3d8-4s2', 'Cu': '3d8-4s2-4p0', 'Zn': '3d10-4s2-4p0', 'Ga': '3d10-4s2-p1', 'Ge': '4s2-4p2', 'As': '3d10-4s2-4p3', 'Se': '4s2-4p4', 'Br': '4s2-4p5', 'Kr': '4s2-4p6', 'Rb': '4s2-4p6-5s1', 'Sr': '4s2-4p6-5s2', 'Y': '4s2-4p6-5s2-4d1', 'Zr': '4s2-4p6-5s2-4d2', 'Nb': '4s2-4p6-5s1-4d4', 'Mo': '4s2-4p6-5s1-4d5', 'Tc': '4s2-4p6-5s2-4d5', 'Ru': '4p6-5s1-4d7', 'Rh': '4p6-5s1-4d8', 'Pd': '4d10-5s0', 'Ag': '4d10-5s1', 'Cd': '4d10-5s2', 'In': '4d10-5s2-5p1', 'Sn': '4d10-5s2-5p2', 'Sb': '5s2-5p3', 'Te': '5s2-5p4', 'I': '5s2-5p5', 'Xe': '5s2-5p6', 'Cs': '5s2-5p6-6s1', 'Ba': '5s2-5p6-6s2', 'La': '5s2-5p6-5d1-6s2', 'Ce': '4f1-5s2-5p6-5d1-6s2', 'Pr': '4f2-5s2-5p6-5d1-6s2', 'Nd': '4f4-5s2-5p6-6s2', 'Pm': '4f5-5s2-5p6-6s2', 'Sm': '4f6-5s2-5p6-6s2', 'Eu': '4f7-5s2-5p6-6s2', 'Gd': '4f7-5s2-5p6-5d1-6s2', 'Tb': '4f9-5s2-5p6-6s2', 'Dy': '4f10-5s2-5p6-6s2', 'Ho': '4f11-5s2-5p6-6s2', 'Er': '4f12-5s2-5p6-6s2', 'Tm': '4f13-5s2-5p6-6s2', 'Yb': '4f14-5s2-5p6-6s2', 'Lu': '4f14-5s2-5p6-5d1-6s2', 'Hf': '5p6-6s2-5d2', 'Ta': '5p6-6s2-5d3', 'W': '5p6-6s2-5d4', 'Re': '5d5-6s2', 'Os': '5d6-6s2', 'Ir': '5d7-6s2', 'Pt': '5d9-6s1', 'Au': '5d10-6s1', 'Hg': '5d10-6s2', 'Tl': '5d10-6s2-6p1', 'Pb': '5d10-6s2-6p2', 'Bi': '5d10-6s2-6p3', 'Po': '5d10-6s2-6p4', 'At': '5d10-6s2-6p5', 'Rn': '6s2-6p6', 'Fr': '6s2-6p6-7s1', 'Ra': '6s2-6p6-7s2', 'Ac': '6s2-6p6-6d1-7s2', 'Th': '6s2-6p6-6d2-7s2', 'Pa': '5f2-6s2-6p6-6d1-7s2', 'U': '5f3-6s2-6p6-6d1-7s2', 'Np': '5f4-6s2-6p6-6d1-7s2', 'Pu': '5f6-6s2-6p6-7s2', 'Am': '5f7-6s2-6p6-7s2', 'Cm': '5f7-6s2-6p6-6d1-7s2', 'Bk': '5f9-6s2-6p6-7s2', 'Cf': '5f10-6s2-6p6-7s2', 'Es': '5f11-6s2-6p6-7s2', 'Fm': '5f12-6s2-6p6-7s2', 'Md': '5f13-6s2-6p6-7s2', 'No': '5f14-6s2-6p6-7s2', 'Lr': '5f14-6s2-6p6-6d1-7s2', 'Rf': '6p6-6d2-7s2', 'Db': '6p6-6d3-7s2', 'Sg': '6p6-6d4-7s2', 'Bh': '6p6-6d5-7s2', 'Hs': '6p6-6d6-7s2', 'Mt': '6p6-6d7-7s2', 'Ds': '6p6-6d8-7s2', 'Rg': '6p6-6d9-7s2', 'Cn': '6p6-6d9-7s2', 'Nh': None, 'Fl': None, 'Mc': None, 'Lv': None, 'Ts': None, 'Og': None}
  periodic_table_dict['num_valence_elec']={'H': 1, 'He': 2, 'Li': 3, 'Be': 2, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8, 'Na': 7, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8, 'K': 9, 'Ca': 10, 'Sc': 11, 'Ti': 12, 'V': 11, 'Cr': 12, 'Mn': 13, 'Fe': 8, 'Co': 9, 'Ni': 10, 'Cu': 11, 'Zn': 12, 'Ga': 13, 'Ge': 4, 'As': 15, 'Se': 6, 'Br': 7, 'Kr': 8, 'Rb': 9, 'Sr': 10, 'Y': 11, 'Zr': 12, 'Nb': 13, 'Mo': 14, 'Tc': 15, 'Ru': 14, 'Rh': 15, 'Pd': 10, 'Ag': 11, 'Cd': 12, 'In': 13, 'Sn': 14, 'Sb': 5, 'Te': 6, 'I': 7, 'Xe': 8, 'Cs': 9, 'Ba': 10, 'La': 11, 'Ce': 12, 'Pr': 13, 'Nd': 14, 'Pm': 15, 'Sm': 16, 'Eu': 17, 'Gd': 18, 'Tb': 19, 'Dy': 20, 'Ho': 21, 'Er': 22, 'Tm': 23, 'Yb': 24, 'Lu': 25, 'Hf': 10, 'Ta': 11, 'W': 12, 'Re': 7, 'Os': 8, 'Ir': 9, 'Pt': 10, 'Au': 11, 'Hg': 12, 'Tl': 13, 'Pb': 14, 'Bi': 15, 'Po': 16, 'At': 17, 'Rn': 8, 'Fr': 9, 'Ra': 10, 'Ac': 11, 'Th': 12, 'Pa': 13, 'U': 14, 'Np': 15, 'Pu': 16, 'Am': 17, 'Cm': 18, 'Bk': 19, 'Cf': 20, 'Es': 21, 'Fm': 22, 'Md': 23, 'No': 24, 'Lr': 25, 'Rf': 10, 'Db': 11, 'Sg': 12, 'Bh': 13, 'Hs': 14, 'Mt': 15, 'Ds': 16, 'Rg': 17, 'Cn': 18, 'Nh': None, 'Fl': None, 'Mc': None, 'Lv': None, 'Ts': None, 'Og': None}
  return periodic_table_dict