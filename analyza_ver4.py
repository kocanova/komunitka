import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from geopy.geocoders import Nominatim
from geopy.distance import geodesic
import time
import os
import json
from datetime import datetime

class KomunitniEnergetika:
    def __init__(self, typ_skupiny="aktivni_zakaznik", verbose=True):
        """
        Inicializace třídy pro optimalizaci komunitní energetiky.
        
        Args:
            typ_skupiny (str): Typ skupiny sdílení ("aktivni_zakaznik" nebo "energeticke_spolecenstvi")
            verbose (bool): Přepínač pro detailní výpisy během běhu programu
        """
        self.verbose = verbose
        self.typ_skupiny = typ_skupiny
        if typ_skupiny not in ["aktivni_zakaznik", "energeticke_spolecenstvi"]:
            raise ValueError("Neplatný typ skupiny sdílení. Povolené hodnoty: 'aktivni_zakaznik', 'energeticke_spolecenstvi'")
        
        self.zdroje_df = None
        self.om_df = None
        self.ceny_df = None
        self.distance_matrix = None
        self.alokacni_klice = None  # Procentuální rozdělení elektřiny
        self.alokovana_energie = None  # Skutečně alokovaná energie v kWh (při použití 15min dat)
        self.priority = None
        self.missing_smart_meters = []
        self.geolocator = Nominatim(user_agent="energy_community_optimizer")
        
        # Nastavení maximálního počtu EAN kódů a teritoriálního omezení podle typu skupiny
        if typ_skupiny == "aktivni_zakaznik":
            self.max_ean = 11
        else:  # energeticke_spolecenstvi
            self.max_ean = 1000
        
        # Nastavení metody alokace - bude určeno podle počtu EAN kódů ve skupině
        self.metoda_alokace = None
        self.max_iterace = None
        
        if self.verbose:
            print(f"Inicializace systému pro optimalizaci komunitní energetiky v ČR (typ: {typ_skupiny})")

    def nacti_data(self, cesta_zdroje, cesta_om, cesta_spotreby_15min=None, cesta_vyroby_15min=None, cesta_ceny=None):
        """
        Načtení dat ze vstupních souborů.
        """
        try:
            self.zdroje_df = pd.read_excel(cesta_zdroje)
            # Převod EAN na řetězec
            self.zdroje_df['EAN'] = self.zdroje_df['EAN'].astype(str)
            if self.verbose:
                print(f"Načteno {len(self.zdroje_df)} zdrojů energie")
                
            self.om_df = pd.read_excel(cesta_om)
            self.om_df['EAN'] = self.om_df['EAN'].astype(str)
            if self.verbose:
                print(f"Načteno {len(self.om_df)} odběrných míst")
                
            # Načtení 15min dat
            self.data_spotreby_15min = None
            self.data_vyroby_15min = None
            if cesta_spotreby_15min:
                self.data_spotreby_15min = pd.read_excel(cesta_spotreby_15min, parse_dates=['Timestamp'], index_col='Timestamp')
                if self.verbose:
                    print(f"Načtena 15min data spotřeby: {self.data_spotreby_15min.shape[0]} časových bodů, {self.data_spotreby_15min.shape[1]} odběrných míst")
            if cesta_vyroby_15min:
                self.data_vyroby_15min = pd.read_excel(cesta_vyroby_15min, parse_dates=['Timestamp'], index_col='Timestamp')
                if self.verbose:
                    print(f"Načtena 15min data výroby: {self.data_vyroby_15min.shape[0]} časových bodů, {self.data_vyroby_15min.shape[1]} zdrojů")
            if cesta_ceny:
                self.ceny_df = pd.read_excel(cesta_ceny)
                # Převod EAN na řetězec pro konzistenci
                self.ceny_df['EAN'] = self.ceny_df['EAN'].astype(str)
                # Kontrola povinných sloupců v cenovém souboru
                required_columns_ceny = ['EAN', 'Silova_cena_kWh', 'Distribuce_cena_kWh']
                missing_columns = [col for col in required_columns_ceny if col not in self.ceny_df.columns]
                if missing_columns:
                    print(f"UPOZORNĚNÍ: Chybí sloupce {missing_columns} v souboru cen. Ekonomická analýza bude méně přesná.")
                if self.verbose:
                    print(f"Načteny cenové údaje pro {len(self.ceny_df)} EAN kódů")
            
            # Kontrola povinných sloupců
            required_columns_zdroje = ['EAN', 'Nazev', 'Adresa', 'Max_vykon', 'Vlastnik']
            required_columns_om = ['EAN', 'Nazev', 'Adresa', 'Max_spotreba', 'Vlastnik', 'Jednotkova_cena']
            
            for col in required_columns_zdroje:
                if col not in self.zdroje_df.columns:
                    raise ValueError(f"Chybí povinný sloupec '{col}' v souboru zdrojů")
            for col in required_columns_om:
                if col not in self.om_df.columns:
                    raise ValueError(f"Chybí povinný sloupec '{col}' v souboru odběrných míst")
            
            # Nastavení metody alokace podle počtu EAN kódů
            self.nastav_metodu_alokace()
            
            # Kontrola limitů pro skupinu sdílení
            if not self.zkontroluj_limity_skupiny():
                print("VAROVÁNÍ: Skupina sdílení nesplňuje limity. Některé funkce mohou být omezeny.")
            
            # Kontrola územních omezení pro energetická společenství
            if self.typ_skupiny == "energeticke_spolecenstvi":
                if not self.zkontroluj_uzemni_omezeni():
                    print("VAROVÁNÍ: Energetické společenství nesplňuje územní omezení (max. 3 ORP).")
            
            return True
                
        except Exception as e:
            print(f"Chyba při načítání dat: {e}")
            return False

    def nastav_metodu_alokace(self):
        """
        Nastaví metodu alokace elektřiny podle velikosti skupiny sdílení.
        """
        total_ean = len(self.zdroje_df) + len(self.om_df) if self.zdroje_df is not None and self.om_df is not None else 0
        
        if total_ean <= 50:
            self.metoda_alokace = "vicekolova_staticka"
            self.max_iterace = 5
        else:
            self.metoda_alokace = "jednokolova_staticka"
            self.max_iterace = 1
        
        if self.verbose:
            print(f"Pro skupinu s {total_ean} EAN kódy nastavena metoda alokace: {self.metoda_alokace}")
            print(f"Maximální počet iterací: {self.max_iterace}")

    def zkontroluj_limity_skupiny(self):
        """
        Kontrola limitů skupiny sdílení podle aktuální legislativy.
        """
        if self.zdroje_df is None or self.om_df is None:
            return True  # Nemůžeme kontrolovat, žádná data nejsou načtena
        
        total_ean = len(self.zdroje_df) + len(self.om_df)
        if self.verbose:
            print(f"Celkový počet EAN kódů ve skupině: {total_ean}")
        
        if self.typ_skupiny == "aktivni_zakaznik":
            if total_ean > 11:
                print(f"UPOZORNĚNÍ: Celkový počet EAN kódů ({total_ean}) přesahuje limit 11 pro aktivního zákazníka")
                return False
        elif self.typ_skupiny == "energeticke_spolecenstvi":
            if total_ean > 1000:
                print(f"UPOZORNĚNÍ: Celkový počet EAN kódů ({total_ean}) přesahuje limit 1000 pro energetické společenství")
                return False
        
        return True

    def zkontroluj_uzemni_omezeni(self):
        """
        Kontrola územních omezení pro energetická společenství.
        """
        if self.typ_skupiny != "energeticke_spolecenstvi" or self.zdroje_df is None or self.om_df is None:
            return True  # Pro aktivní zákazníky neplatí územní omezení nebo nemáme data
        
        # Získání unikátních adres
        adresy = set(self.zdroje_df['Adresa'].tolist() + self.om_df['Adresa'].tolist())
        
        # Extrakce měst z adres (zjednodušený model - předpokládáme, že města jsou oddělena čárkou)
        mesta = set()
        for adresa in adresy:
            if isinstance(adresa, str) and ',' in adresa:
                mesto = adresa.split(',')[-1].strip()
                if mesto:
                    mesta.add(mesto)
        
        # Zjednodušená implementace - v reálné aplikaci by bylo potřeba použít databázi ORP
        if len(mesta) > 3:
            if self.verbose:
                print(f"UPOZORNĚNÍ: Účastníci skupiny sdílení se nachází ve více než 3 městech: {mesta}")
            return False
        
        if self.verbose:
            print(f"Účastníci skupiny sdílení se nachází ve městech: {mesta}")
        
        return True

    def zkontroluj_chytre_merice(self, smart_meter_sloupec='Chytry_meric'):
        """
        Kontrola přítomnosti chytrých měřičů u zdrojů a odběrných míst.
        """
        self.missing_smart_meters = []
        if smart_meter_sloupec not in self.zdroje_df.columns:
            self.zdroje_df[smart_meter_sloupec] = True
        if smart_meter_sloupec not in self.om_df.columns:
            self.om_df[smart_meter_sloupec] = True
        for idx, row in self.zdroje_df.iterrows():
            if not row[smart_meter_sloupec]:
                self.missing_smart_meters.append({
                    'EAN': row['EAN'],
                    'Typ': 'Zdroj',
                    'Nazev': row['Nazev'],
                    'Adresa': row['Adresa']
                })
        for idx, row in self.om_df.iterrows():
            if not row[smart_meter_sloupec]:
                self.missing_smart_meters.append({
                    'EAN': row['EAN'],
                    'Typ': 'Odběrné místo',
                    'Nazev': row['Nazev'],
                    'Adresa': row['Adresa']
                })
        if self.verbose:
            print(f"Kontrola chytrých měřičů dokončena. Nalezeno {len(self.missing_smart_meters)} míst bez chytrého měřiče.")
def vypocti_vzdalenosti(self):
    """
    Výpočet vzdáleností mezi zdroji a odběrnými místy.
    """
    if self.zdroje_df is None or self.om_df is None:
        print("Nejsou načtena potřebná data")
        return
    self.distance_matrix = pd.DataFrame(
        index=self.om_df['EAN'],
        columns=self.zdroje_df['EAN']
    )
    if self.verbose:
        print("Získávání souřadnic pro adresy (může trvat delší dobu)...")
    zdroje_coords = {}
    om_coords = {}
    for idx, row in self.zdroje_df.iterrows():
        coords = self._get_coordinates(row['Adresa'])
        zdroje_coords[row['EAN']] = coords
        time.sleep(1)
    for idx, row in self.om_df.iterrows():
        coords = self._get_coordinates(row['Adresa'])
        om_coords[row['EAN']] = coords
        time.sleep(1)
    for om_ean, om_coord in om_coords.items():
        for zdroj_ean, zdroj_coord in zdroje_coords.items():
            if om_coord and zdroj_coord:
                distance = geodesic(om_coord, zdroj_coord).kilometers
                self.distance_matrix.loc[om_ean, zdroj_ean] = round(distance, 2)
            else:
                self.distance_matrix.loc[om_ean, zdroj_ean] = np.nan
    if self.verbose:
        print("Výpočet vzdáleností dokončen.")
    def jsou_ve_stejne_lokaci(self, ean1, ean2):
        """
        Kontroluje, zda jsou dva EAN kódy ve stejné lokaci (za stejnou hlavní domovní skříní).
        
        V této zjednodušené verzi kontrolujeme pouze shodu adres, ale v reálném nasazení
        by bylo potřeba přesnější kontroly, např. identifikace hlavní domovní skříně.
        """
        # Najít adresy pro oba EAN kódy
        adresa1 = None
        adresa2 = None
        
        # Hledat v zdrojích
        for _, zdroj_row in self.zdroje_df.iterrows():
            if zdroj_row['EAN'] == ean1:
                adresa1 = zdroj_row['Adresa']
            if zdroj_row['EAN'] == ean2:
                adresa2 = zdroj_row['Adresa']
        
        # Hledat v odběrných místech
        for _, om_row in self.om_df.iterrows():
            if om_row['EAN'] == ean1:
                adresa1 = om_row['Adresa']
            if om_row['EAN'] == ean2:
                adresa2 = om_row['Adresa']
        
        # Kontrola shody adres (v reálném nasazení by bylo potřeba komplexnější kontroly)
        if adresa1 is not None and adresa2 is not None:
            return adresa1 == adresa2
        
        return False

    def vytvor_priority(self):
        """
        Vytvoří matici priorit pro každé odběrné místo podle vzdálenosti, adresy a vlastníka.
        Kontroluje také limit maximálně 5 výroben sdílejících do jednoho OM.
        """
        if self.distance_matrix is None:
            print("Nejsou spočítány vzdálenosti")
            return
        self.priority = pd.DataFrame(
            index=self.om_df['EAN'],
            columns=self.zdroje_df['EAN']
        )
        for om_idx, om_row in self.om_df.iterrows():
            om_ean = om_row['EAN']
            om_adresa = om_row['Adresa']
            om_vlastnik = om_row['Vlastnik']
            zdroj_metriky = []
            for zdroj_idx, zdroj_row in self.zdroje_df.iterrows():
                zdroj_ean = zdroj_row['EAN']
                zdroj_adresa = zdroj_row['Adresa']
                zdroj_vlastnik = zdroj_row['Vlastnik']
                distance = self.distance_matrix.loc[om_ean, zdroj_ean]
                same_address = om_adresa == zdroj_adresa
                same_owner = om_vlastnik == zdroj_vlastnik
                priority_score = distance
                if same_address:
                    priority_score = -1000
                elif same_owner:
                    priority_score = -500 + distance
                zdroj_metriky.append({
                    'zdroj_ean': zdroj_ean,
                    'distance': distance,
                    'same_address': same_address,
                    'same_owner': same_owner,
                    'priority_score': priority_score
                })
            sorted_zdroje = sorted(zdroj_metriky, key=lambda x: x['priority_score'])
            for idx, zdroj_info in enumerate(sorted_zdroje):
                self.priority.loc[om_ean, zdroj_info['zdroj_ean']] = idx + 1
        
        # Kontrola limitu počtu výroben na jedno OM
        for om_ean in self.om_df['EAN']:
            prioritni_zdroje = self.priority.loc[om_ean].sort_values()
            if len(prioritni_zdroje) > 5:
                if self.verbose:
                    print(f"UPOZORNĚNÍ: OM {om_ean} má přiřazeno {len(prioritni_zdroje)} zdrojů, což přesahuje limit 5 výroben na jedno OM")
                # Ponecháme pouze 5 zdrojů s nejvyšší prioritou (nejnižší hodnotou)
                zdroje_k_odstraneni = prioritni_zdroje.iloc[5:].index
                for zdroj_ean in zdroje_k_odstraneni:
                    self.priority.loc[om_ean, zdroj_ean] = float('inf')  # Nastavit nekonečnou prioritu = nebude používán
                if self.verbose:
                    print(f"Pro OM {om_ean} ponecháno pouze 5 zdrojů s nejvyšší prioritou")
        
        if self.verbose:
            print("Vytvoření priorit dokončeno.")
    def optimalizuj_alokacni_klice(self, iterace=None, pouzit_15min_data=False):
        """
        Optimalizace alokačních klíčů pro maximalizaci celkových úspor.
        """
        if self.priority is None:
            print("Nejsou vytvořeny priority")
            return
        
        if self.verbose:
            print("Spouštím optimalizaci alokačních klíčů pro maximalizaci úspor...")
        
        # Inicializace alokačních klíčů
        self.alokacni_klice = pd.DataFrame(0.0, index=self.zdroje_df['EAN'], columns=self.om_df['EAN'])
        
        # Pro každý zdroj vytvoříme seznam prioritizovaných OM podle potenciální úspory
        for zdroj_idx, zdroj_row in self.zdroje_df.iterrows():
            zdroj_ean = zdroj_row['EAN']
            
            # Seznam potenciálních příjemců s informacemi o potenciální úspoře
            om_info = []
            
            for om_idx, om_row in self.om_df.iterrows():
                om_ean = om_row['EAN']
                
                # Získání cenových údajů
                ceny_om = self.ceny_df.loc[self.ceny_df['EAN'] == om_ean] if self.ceny_df is not None else None
                if len(ceny_om) > 0:
                    cena_silova = ceny_om.iloc[0]['Silova_cena_kWh']
                    cena_distribuce = ceny_om.iloc[0]['Distribuce_cena_kWh']
                else:
                    cena_silova = om_row['Jednotkova_cena'] * 0.6
                    cena_distribuce = om_row['Jednotkova_cena'] * 0.4
                
                # Zjistíme, zda jsou zdroj a OM ve stejné lokaci
                same_location = self.jsou_ve_stejne_lokaci(zdroj_ean, om_ean)
                
                # Vypočteme potenciální úsporu na kWh
                potencialni_uspora = cena_silova
                if same_location:
                    potencialni_uspora += cena_distribuce
                
                # Zohledníme také prioritu (vzdálenost, vlastnictví) s nižší váhou
                priorita = self.priority.loc[om_ean, zdroj_ean] if om_ean in self.priority.index and zdroj_ean in self.priority.columns else float('inf')
                if priorita != float('inf'):
                    # Přidáme bonus za vysokou prioritu (max 20% z ceny)
                    priorita_bonus = cena_silova * 0.2 * (10 - min(priorita, 9)) / 9
                    potencialni_uspora += priorita_bonus
                
                # Omezení na 5 výroben na OM
                if priorita <= 5:
                    om_info.append({
                        'om_ean': om_ean,
                        'uspora_per_kwh': potencialni_uspora,
                        'max_spotreba': om_row['Max_spotreba'],
                        'priorita': priorita
                    })
            
            # Seřadíme odběrná místa podle potenciální úspory (sestupně)
            om_info = sorted(om_info, key=lambda x: x['uspora_per_kwh'], reverse=True)
            
            # Postupně alokujeme energii, začínáme od OM s nejvyšší potenciální úsporou
            remaining_capacity = 1.0  # 100% kapacity zdroje
            for om in om_info:
                if remaining_capacity <= 0:
                    break
                    
                # Určíme alokační procento (omezeno zbývající kapacitou)
                # Pro jednoduchost rozdělujeme kapacitu rovnoměrně mezi horních N odběrných míst
                # V reálné aplikaci by bylo vhodné sofistikovanější přiřazení podle množství spotřeby
                alloc_percent = min(0.2, remaining_capacity)  # Max 20% na jedno OM, aby bylo více diverzifikováno
                
                self.alokacni_klice.loc[zdroj_ean, om['om_ean']] = alloc_percent
                remaining_capacity -= alloc_percent
        
        # Pokud máme 15minutová data, použijeme je pro simulaci skutečných toků energie
        if pouzit_15min_data and self.data_spotreby_15min is not None and self.data_vyroby_15min is not None:
            # Simulace alokace v 15min intervalech (pro ekonomickou analýzu)
            self.simuluj_alokaci_s_15min_daty(iterace if iterace is not None else self.max_iterace)
        else:
            # Pro případy bez 15min dat vytvoříme aproximaci alokované energie
            self.simuluj_alokaci_bez_15min_dat()
            
        if self.verbose:
            print("Optimalizace alokačních klíčů dokončena.")
    def simuluj_alokaci_s_15min_daty(self, iterace):
        """
        Simuluje alokaci elektřiny s využitím 15‑minutových dat pro ekonomickou analýzu.
        Používá již nastavené alokační klíče a priority.
        """
        if self.verbose:
            print("Simulace alokace s využitím 15‑minutových dat...")
        
        # Inicializujeme matici pro skutečně alokovanou energii
        self.alokovana_energie = pd.DataFrame(0.0, index=self.zdroje_df['EAN'], columns=self.om_df['EAN'])
        
        total_intervals = len(self.data_spotreby_15min)
        print(f"Simuluji {total_intervals} 15‑minutových intervalů...")
        
        for timestamp in self.data_spotreby_15min.index:
            print(f"Interval: {timestamp}")
            current_consumption = self.data_spotreby_15min.loc[timestamp].copy()
            current_production = self.data_vyroby_15min.loc[timestamp].copy()
            
            # Pro každý interval provést iterace alokací podle alokačních klíčů
            remaining_production = current_production.copy()
            remaining_consumption = current_consumption.copy()
            
            for it in range(iterace):
                print(f"  Iterace {it+1}/{iterace} pro interval {timestamp}")
                
                # Temporary storage for this iteration's allocation
                iteration_allocation = pd.DataFrame(0.0, index=self.zdroje_df['EAN'], columns=self.om_df['EAN'])
                
                # First pass: allocate according to percentages
                for zdroj_ean in self.zdroje_df['EAN']:
                    if zdroj_ean not in remaining_production.index or remaining_production[zdroj_ean] <= 0:
                        continue
                    
                    for om_ean in self.om_df['EAN']:
                        if om_ean not in remaining_consumption.index or remaining_consumption[om_ean] <= 0:
                            continue
                        
                        # Get allocation percentage from the allocation matrix
                        alloc_percent = self.alokacni_klice.loc[zdroj_ean, om_ean]
                        if alloc_percent <= 0:
                            continue
                        
                        # Calculate how much energy to allocate
                        energy_to_allocate = remaining_production[zdroj_ean] * alloc_percent
                        energy_allocated = min(energy_to_allocate, remaining_consumption[om_ean])
                        
                        # Update remaining values
                        remaining_production[zdroj_ean] -= energy_allocated
                        remaining_consumption[om_ean] -= energy_allocated
                        
                        # Record this allocation
                        iteration_allocation.loc[zdroj_ean, om_ean] = energy_allocated
                
                # Add this iteration's allocation to the total
                self.alokovana_energie += iteration_allocation
                
                # If nothing left to allocate or consume, break
                if (remaining_production <= 0).all() or (remaining_consumption <= 0).all():
                    break
        
        if self.verbose:
            print("Simulace alokace s 15‑minutovými daty dokončena.")
            print(f"Celková alokovaná energie: {self.alokovana_energie.sum().sum():.2f} kWh")
    def simuluj_alokaci_bez_15min_dat(self):
        """
        Simuluje alokaci elektřiny bez 15-minutových dat.
        Vytváří aproximaci alokované energie na základě maximálních hodnot.
        """
        if self.verbose:
            print("Simulace alokace bez 15‑minutových dat...")
        
        # Inicializujeme matici pro alokovanou energii
        self.alokovana_energie = pd.DataFrame(0.0, index=self.zdroje_df['EAN'], columns=self.om_df['EAN'])
        
        # Použijeme podobný postup jako v původní metodě optimalizuj_alokacni_klice
        remaining_capacity = self.zdroje_df['Max_vykon'].copy()
        remaining_consumption = self.om_df['Max_spotreba'].copy()
        
        for it in range(self.max_iterace):
            print(f"Iterace {it+1}/{self.max_iterace}")
            
            # Temporary storage for this iteration's allocation
            iteration_allocation = pd.DataFrame(0.0, index=self.zdroje_df['EAN'], columns=self.om_df['EAN'])
            
            for zdroj_idx, zdroj_row in self.zdroje_df.iterrows():
                zdroj_ean = zdroj_row['EAN']
                if remaining_capacity[zdroj_idx] <= 0:
                    continue
                
                for om_idx, om_row in self.om_df.iterrows():
                    om_ean = om_row['EAN']
                    if remaining_consumption[om_idx] <= 0:
                        continue
                    
                    # Get allocation percentage
                    alloc_percent = self.alokacni_klice.loc[zdroj_ean, om_ean]
                    if alloc_percent <= 0:
                        continue
                    
                    # Calculate allocation
                    energy_to_allocate = remaining_capacity[zdroj_idx] * alloc_percent
                    energy_allocated = min(energy_to_allocate, remaining_consumption[om_idx])
                    
                    # Update remaining values
                    remaining_capacity[zdroj_idx] -= energy_allocated
                    remaining_consumption[om_idx] -= energy_allocated
                    
                    # Record this allocation
                    iteration_allocation.loc[zdroj_ean, om_ean] = energy_allocated
            
            # Add this iteration's allocation to the total
            self.alokovana_energie += iteration_allocation
            
            # If nothing left to allocate or consume, break
            if (remaining_capacity <= 0).all() or (remaining_consumption <= 0).all():
                break
        
        if self.verbose:
            print("Simulace alokace bez 15‑minutových dat dokončena.")
            print(f"Celková alokovaná energie: {self.alokovana_energie.sum().sum():.2f} kWh")

    def _get_coordinates(self, address, max_attempts=3):
        """
        Získání geografických souřadnic pro adresu.
        """
        for attempt in range(max_attempts):
            try:
                address = str(address).strip()
                if not address:
                    return None
                location = self.geolocator.geocode(address + ", Czech Republic")
                if location:
                    return (location.latitude, location.longitude)
                else:
                    location = self.geolocator.geocode(address)
                    return (location.latitude, location.longitude) if location else None
            except Exception as e:
                if self.verbose:
                    print(f"Chyba geokódování pro {address}: {e}")
                time.sleep(1)
                continue
        return None

    def vypocti_ekonomickou_analyzu(self):
        """
        Ekonomická analýza úspor při sdílení energie.
        """
        if self.alokacni_klice is None:
            print("Nejsou vytvořeny alokační klíče")
            return None
        
        if self.ceny_df is None:
            print("UPOZORNĚNÍ: Nejsou k dispozici cenové údaje, analýza bude méně přesná")
            return self._vypocti_ekonomickou_analyzu_bez_cen()
        
        # Použijeme skutečně alokovanou energii, pokud je k dispozici
        energie_matice = self.alokovana_energie if hasattr(self, 'alokovana_energie') and self.alokovana_energie is not None else None
        
        if energie_matice is None:
            print("UPOZORNĚNÍ: Není k dispozici matice alokované energie, nelze provést ekonomickou analýzu")
            return None
        
        if self.verbose:
            print("Provádím ekonomickou analýzu s použitím přesných cenových údajů.")
        
        results = []
        for om_idx, om_row in self.om_df.iterrows():
            om_ean = om_row['EAN']
            om_adresa = om_row['Adresa']
            
            # Pokud jsou 15‑minutová data k dispozici, celková spotřeba se počítá ze součtu spotřeby
            if self.data_spotreby_15min is not None and om_ean in self.data_spotreby_15min.columns:
                om_spotreba = self.data_spotreby_15min[om_ean].sum()
            else:
                om_spotreba = om_row['Max_spotreba']
            
            # Získání cenových údajů z cenového souboru
            ceny_om = self.ceny_df.loc[self.ceny_df['EAN'] == om_ean]
            if len(ceny_om) == 0:
                print(f"UPOZORNĚNÍ: Pro OM {om_ean} nejsou k dispozici cenové údaje, používám standardní ceny")
                cena_silova = om_row['Jednotkova_cena'] * 0.6  # Přibližný odhad silové ceny
                cena_distribuce = om_row['Jednotkova_cena'] * 0.4  # Přibližný odhad distribuční ceny
            else:
                # Používáme názvy sloupců z vašeho souboru 'Ceny.xlsx'
                cena_silova = ceny_om.iloc[0]['Silova_cena_kWh']
                cena_distribuce = ceny_om.iloc[0]['Distribuce_cena_kWh']
            
            celkova_energie = 0
            celkova_energie_bez_ds = 0  # Elektřina sdílená bez využití distribuční soustavy
            uspora_silova = 0
            uspora_distribuce = 0
            
            for zdroj_idx, zdroj_row in self.zdroje_df.iterrows():
                zdroj_ean = zdroj_row['EAN']
                zdroj_adresa = zdroj_row['Adresa']
                
                # Použijeme matici alokované energie
                allocated_energy = energie_matice.loc[zdroj_ean, om_ean]
                
                if allocated_energy > 0:
                    celkova_energie += allocated_energy
                    
                    # Kontrola, zda jsou zdroj a OM ve stejném místě (za stejnou HDV)
                    same_location = self.jsou_ve_stejne_lokaci(om_ean, zdroj_ean)
                    
                    # Úspora na silové složce se započítá vždy
                    uspora_silova += allocated_energy * cena_silova
                    
                    if same_location:
                        # Elektřina sdílená bez využití distribuční soustavy
                        celkova_energie_bez_ds += allocated_energy
                        
                        # Úspora na distribuci - používáme přesnou cenu z cenového souboru
                        uspora_distribuce += allocated_energy * cena_distribuce
                        
                        if self.verbose:
                            print(f"OM {om_ean} a zdroj {zdroj_ean} jsou ve stejné lokaci - úspora zahrnuje i distribuční poplatky")
            
            zbyvajici_energie = max(0, om_spotreba - celkova_energie)
            celkova_uspora = uspora_silova + uspora_distribuce
            
            # Výpočet celkové ceny, kterou by zákazník zaplatil bez sdílení
            cena_bez_sdileni = om_spotreba * (cena_silova + cena_distribuce)
            
            results.append({
                'EAN': om_ean,
                'Nazev': om_row['Nazev'],
                'Adresa': om_adresa,
                'Celkova_spotreba': om_spotreba,
                'Sdilena_energie_celkem': celkova_energie,
                'Sdilena_energie_bez_DS': celkova_energie_bez_ds,
                'Zbyvajici_energie': zbyvajici_energie,
                'Uspora_silova': uspora_silova,
                'Uspora_distribuce': uspora_distribuce,
                'Celkova_uspora': celkova_uspora,
                'Procento_uspory': (celkova_uspora / cena_bez_sdileni * 100) if cena_bez_sdileni > 0 else 0,
                'Cena_silova': cena_silova,
                'Cena_distribuce': cena_distribuce
            })
        
        return pd.DataFrame(results)

    def _vypocti_ekonomickou_analyzu_bez_cen(self):
        """
        Záložní metoda pro výpočet ekonomické analýzy, pokud nejsou k dispozici cenové údaje.
        Používá přibližný odhad podílu silové a distribuční složky v ceně elektřiny.
        """
        if self.verbose:
            print("Struktura ceny elektřiny (přibližné hodnoty):")
            print("- Silová elektřina (ušetříte při sdílení): ~27% z celkové ceny")
            print("- Distribuce (ušetříte při sdílení v rámci stejné lokace): ~39% z celkové ceny")
            print("- Příspěvek na podporované zdroje: ~7% z celkové ceny")
            print("- Ostatní poplatky a daně: ~27% z celkové ceny")
        
        # Použijeme skutečně alokovanou energii, pokud je k dispozici
        energie_matice = self.alokovana_energie if hasattr(self, 'alokovana_energie') and self.alokovana_energie is not None else None
        
        if energie_matice is None:
            print("UPOZORNĚNÍ: Není k dispozici matice alokované energie, nelze provést ekonomickou analýzu")
            return None
        
        results = []
        for om_idx, om_row in self.om_df.iterrows():
            om_ean = om_row['EAN']
            om_adresa = om_row['Adresa']
            
            # Pokud jsou 15‑minutová data k dispozici, celková spotřeba se počítá ze součtu spotřeby
            if self.data_spotreby_15min is not None and om_ean in self.data_spotreby_15min.columns:
                om_spotreba = self.data_spotreby_15min[om_ean].sum()
            else:
                om_spotreba = om_row['Max_spotreba']
            
            om_cena = om_row['Jednotkova_cena']
            cena_silova = om_cena * 0.6  # Přibližný odhad silové ceny
            cena_distribuce = om_cena * 0.4  # Přibližný odhad distribuční ceny
            
            celkova_energie = 0
            celkova_energie_bez_ds = 0
            uspora_silova = 0
            uspora_distribuce = 0
            
            for zdroj_idx, zdroj_row in self.zdroje_df.iterrows():
                zdroj_ean = zdroj_row['EAN']
                zdroj_adresa = zdroj_row['Adresa']
                
                # Použijeme matici alokované energie
                allocated_energy = energie_matice.loc[zdroj_ean, om_ean]
                
                if allocated_energy > 0:
                    celkova_energie += allocated_energy
                    
                    # Kontrola, zda jsou zdroj a OM ve stejném místě
                    same_location = self.jsou_ve_stejne_lokaci(om_ean, zdroj_ean)
                    
                    # Úspora na silové složce se započítá vždy
                    uspora_silova += allocated_energy * cena_silova
                    
                    if same_location:
                        celkova_energie_bez_ds += allocated_energy
                        uspora_distribuce += allocated_energy * cena_distribuce
            
            zbyvajici_energie = max(0, om_spotreba - celkova_energie)
            celkova_uspora = uspora_silova + uspora_distribuce
            
            # Výpočet celkové ceny, kterou by zákazník zaplatil bez sdílení
            cena_bez_sdileni = om_spotreba * (cena_silova + cena_distribuce)
            
            results.append({
                'EAN': om_ean,
                'Nazev': om_row['Nazev'],
                'Adresa': om_adresa,
                'Celkova_spotreba': om_spotreba,
                'Sdilena_energie_celkem': celkova_energie,
                'Sdilena_energie_bez_DS': celkova_energie_bez_ds,
                'Zbyvajici_energie': zbyvajici_energie,
                'Uspora_silova': uspora_silova,
                'Uspora_distribuce': uspora_distribuce,
                'Celkova_uspora': celkova_uspora,
                'Procento_uspory': (celkova_uspora / cena_bez_sdileni * 100) if cena_bez_sdileni > 0 else 0
            })
        
        return pd.DataFrame(results)

    def generuj_data_pro_edc(self, output_dir='vystupy'):
        """
        Generuje data pro registraci u Elektroenergetického datového centra (EDC).
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Vygenerovat data pro registraci zdrojů (EANd)
        eand_data = []
        for _, zdroj_row in self.zdroje_df.iterrows():
            eand_data.append({
                'EAN': zdroj_row['EAN'],
                'Nazev': zdroj_row['Nazev'],
                'Adresa': zdroj_row['Adresa'],
                'Max_vykon': float(zdroj_row['Max_vykon']),
                'Vlastnik': zdroj_row['Vlastnik']
            })
        
        # Vygenerovat data pro registraci odběrných míst (EANo)
        eano_data = []
        for _, om_row in self.om_df.iterrows():
            om_ean = om_row['EAN']
            om_adresa = om_row['Adresa']
            
            if self.priority is None:
                print("UPOZORNĚNÍ: Nejsou vytvořeny priority, nelze generovat kompletní data pro EDC")
                zdroje_info = []
            else:
                prioritni_zdroje = self.priority.loc[om_ean].sort_values()
                # Omezit na 5 nejvyšších priorit
                prioritni_zdroje = prioritni_zdroje.iloc[:5] if len(prioritni_zdroje) > 5 else prioritni_zdroje
                
                zdroje_info = []
                for zdroj_ean, priorita in prioritni_zdroje.items():
                    if priorita != float('inf'):  # Pouze platné priority
                        # Alokační procento je přímo z matice alokačních klíčů (v procentech)
                        alokace_procento = self.alokacni_klice.loc[zdroj_ean, om_ean] * 100
                        
                        # Kontrola, zda jsou zdroj a OM ve stejném místě
                        same_location = self.jsou_ve_stejne_lokaci(om_ean, zdroj_ean)
                        
                        zdroje_info.append({
                            'EANd': zdroj_ean,
                            'Priorita': int(priorita),
                            'Alokace_procento': float(alokace_procento),
                            'Vyuziti_DS': not same_location  # False pokud jsou ve stejném místě, True jinak
                        })
            
            eano_data.append({
                'EANo': om_ean,
                'Nazev': om_row['Nazev'],
                'Adresa': om_adresa,
                'Zdroje': zdroje_info
            })
        
        # Uložit data do JSON souboru
        edc_data = {
            'Typ_skupiny': self.typ_skupiny,
            'EANd': eand_data,
            'EANo': eano_data
        }
        
        filename = os.path.join(output_dir, f"{timestamp}_data_pro_edc.json")
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(edc_data, f, ensure_ascii=False, indent=2)
        
        if self.verbose:
            print(f"Data pro registraci u EDC uložena do: {filename}")

    def generuj_vystupy(self, output_dir='vystupy'):
        """
        Generování výstupů modelu.
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if self.distance_matrix is not None:
            filename = os.path.join(output_dir, f"{timestamp}_vzdalenosti.xlsx")
            self.distance_matrix.to_excel(filename)
            if self.verbose:
                print(f"Matice vzdáleností uložena do: {filename}")
        
        if self.priority is not None:
            filename = os.path.join(output_dir, f"{timestamp}_priority.xlsx")
            self.priority.to_excel(filename)
            if self.verbose:
                print(f"Matice priorit uložena do: {filename}")
        
        if self.alokacni_klice is not None:
            filename = os.path.join(output_dir, f"{timestamp}_alokacni_klice.xlsx")
            self.alokacni_klice.to_excel(filename)
            if self.verbose:
                print(f"Alokační klíče uloženy do: {filename}")
        
        if hasattr(self, 'alokovana_energie') and self.alokovana_energie is not None:
            filename = os.path.join(output_dir, f"{timestamp}_alokovana_energie.xlsx")
            self.alokovana_energie.to_excel(filename)
            if self.verbose:
                print(f"Alokovaná energie uložena do: {filename}")
        
        ekonomicka_analyza = self.vypocti_ekonomickou_analyzu()
        if ekonomicka_analyza is not None:
            filename = os.path.join(output_dir, f"{timestamp}_ekonomicka_analyza.xlsx")
            ekonomicka_analyza.to_excel(filename, index=False)
            if self.verbose:
                print(f"Ekonomická analýza uložena do: {filename}")
        
        if self.missing_smart_meters:
            filename = os.path.join(output_dir, f"{timestamp}_chybejici_chytre_merice.xlsx")
            pd.DataFrame(self.missing_smart_meters).to_excel(filename, index=False)
            if self.verbose:
                print(f"Seznam chybějících chytrých měřičů uložen do: {filename}")
        
        self._vytvor_vizualizaci(os.path.join(output_dir, f"{timestamp}_vizualizace.png"))
        
        # Vygenerovat data pro EDC
        self.generuj_data_pro_edc(output_dir)

    def _vytvor_vizualizaci(self, filename):
        """
        Vytvoření vizualizace toků energie v komunitě.
        """
        if self.alokacni_klice is None:
            return
        
        plt.figure(figsize=(12, 8))
        plt.title("Toky energie v komunitě")
        plt.xlabel("Odběrná místa")
        plt.ylabel("Zdroje")
        im = plt.imshow(self.alokacni_klice.values, cmap='Blues')
        plt.xticks(range(len(self.alokacni_klice.columns)), self.alokacni_klice.columns, rotation=90)
        plt.yticks(range(len(self.alokacni_klice.index)), self.alokacni_klice.index)
        plt.colorbar(im, label="Alokační klíč")
        plt.tight_layout()
        plt.savefig(filename)
        if self.verbose:
            print(f"Vizualizace uložena do: {filename}")


def main():
    print("\n===== OPTIMALIZACE KOMUNITNÍ ENERGETIKY =====\n")
    
    print("Vyberte typ skupiny sdílení:")
    print("1. Aktivní zákazník (max. 11 EAN)")
    print("2. Energetické společenství (max. 1000 EAN)")
    volba = input("Váš výběr (1/2): ")
    
    typ_skupiny = "aktivni_zakaznik" if volba == "1" else "energeticke_spolecenstvi"
    
    optimizer = KomunitniEnergetika(typ_skupiny=typ_skupiny, verbose=True)
    
    print("\nPro spuštění optimalizace potřebuji následující vstupní soubory:")
    print("1. Excel soubor se zdroji energie (povinný)")
    print("2. Excel soubor s odběrnými místy (povinný)")
    print("3. Excel soubor s 15‑minutovými daty spotřeby (doporučeno)")
    print("4. Excel soubor s 15‑minutovými daty výroby (doporučeno)")
    print("5. Volitelně Excel soubor s cenami energie")
    
    print("\nProsím zadejte cesty k souborům:")
    while True:
        zdroje_path = input("Cesta k souboru se zdroji: ")
        om_path = input("Cesta k souboru s odběrnými místy: ")
        
        # Základní validace vstupních souborů
        if not os.path.exists(zdroje_path):
            print(f"CHYBA: Soubor {zdroje_path} nebyl nalezen. Zkuste znovu.")
            continue
        
        if not os.path.exists(om_path):
            print(f"CHYBA: Soubor {om_path} nebyl nalezen. Zkuste znovu.")
            continue
        
        # Volitelné soubory
        spotreba_15min_path = input("Cesta k Excel souboru s 15‑minutovými daty spotřeby (pokud nemáte, stiskněte Enter): ")
        vyroba_15min_path = input("Cesta k Excel souboru s 15‑minutovými daty výroby (pokud nemáte, stiskněte Enter): ")
        ceny_path = input("Cesta k souboru s cenami (pokud nemáte, stiskněte Enter): ")
        
        # Validace volitelných souborů
        if spotreba_15min_path and not os.path.exists(spotreba_15min_path):
            print(f"CHYBA: Soubor {spotreba_15min_path} nebyl nalezen. Zkuste znovu.")
            continue
        
        if vyroba_15min_path and not os.path.exists(vyroba_15min_path):
            print(f"CHYBA: Soubor {vyroba_15min_path} nebyl nalezen. Zkuste znovu.")
            continue
        
        if ceny_path and not os.path.exists(ceny_path):
            print(f"CHYBA: Soubor {ceny_path} nebyl nalezen. Zkuste znovu.")
            continue
        
        # Nastavení None pro prázdné volitelné soubory
        spotreba_15min_path = spotreba_15min_path if spotreba_15min_path else None
        vyroba_15min_path = vyroba_15min_path if vyroba_15min_path else None
        ceny_path = ceny_path if ceny_path else None
        
        # Pokus o načtení dat
        try:
            if optimizer.nacti_data(zdroje_path, om_path, spotreba_15min_path, vyroba_15min_path, ceny_path):
                break
            else:
                print("Nepodařilo se načíst data. Zkuste znovu.")
        except Exception as e:
            print(f"Chyba při načítání dat: {e}")
            print("Zkuste znovu zadat soubory.")
    
    # Příprava optimalizace
    optimizer.zkontroluj_chytre_merice()
    optimizer.vypocti_vzdalenosti()
    optimizer.vytvor_priority()
    
    iterace = optimizer.max_iterace
    
    # Dotaz na použití 15‑minutových dat
    pouzit_15min = False
    if spotreba_15min_path and vyroba_15min_path:
        while True:
            odpoved = input("\nChcete použít 15‑minutová data pro přesnější optimalizaci? (ano/ne): ").lower()
            if odpoved in ['ano', 'a', 'yes', 'y']:
                pouzit_15min = True
                break
            elif odpoved in ['ne', 'n', 'no']:
                pouzit_15min = False
                break
            else:
                print("Neplatná odpověď. Použijte 'ano' nebo 'ne'.")
    
    # Optimalizace alokačních klíčů
    optimizer.optimalizuj_alokacni_klice(iterace, pouzit_15min)
    
    # Volba výstupního adresáře
    while True:
        output_dir = input("\nZadejte adresář pro výstupy (výchozí: 'vystupy'): ")
        if not output_dir:
            output_dir = 'vystupy'
        
        try:
            # Pokus o vytvoření adresáře, pokud neexistuje
            os.makedirs(output_dir, exist_ok=True)
            break
        except Exception as e:
            print(f"Nepodařilo se vytvořit adresář {output_dir}. Chyba: {e}")
            print("Zkuste zadat jiný adresář.")
    
    # Generování výstupů
    optimizer.generuj_vystupy(output_dir)
    
    print("\n===== OPTIMALIZACE DOKONČENA =====")
    print(f"Výsledky byly uloženy do adresáře: {output_dir}")

# Umožnit spuštění skriptu přímo
if __name__ == "__main__":
    main()
           