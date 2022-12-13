import numpy as np
from astropy.io import fits
import inspect
import time

def read_fits_table2dict(f,hdu):
    '''
    reads arbitrary fits table and creates dictionary with information
    using the first column as identifier
    '''
    
    hdulist = fits.open(f)
    tbdata = hdulist[hdu].data
    tbcols = hdulist[hdu].columns
    hdulist.close()
    info = {}
    for col in tbcols.names:
        info[col] = tbdata[col] #create dictionary from table columns

    # find the column containing the IDs
    try:
        if 'ID' in list(info[tbcols.names[0]]):
            ids = info[tbcols.names[0]]
        else:
            try: 
                ids = info['ID']
            except:
                here_ids = [i for i in tbcols.names if 'ID' in i]
                ids = info[here_ids[0]]    
    except:
        ids = info[tbcols.names[0]]

    info_ids = {}
    for i in range(len(ids)):
        id_params = {}
        for c in info.keys():
            id_params[c] = info[c][i]
        info_ids[str(ids[i])] = id_params

    return info,info_ids


def write_table_from_dict_in_dict(d,f):
    '''
    Write a fits table using a dictionary with dictionaries using the
    keys as first column and the keys of the inner dictionaries as other
    columns. 
    d = dictionary that is gonna be saved in the output file
    f = name and directory of output file
    '''

    all_cols_content = reshape_dict_in_dict(d,'ID')
    # test datatype of contents
    data_types = {}
    for c in all_cols_content: #c is each key of the dictionary (i.e. each column name)
        values = all_cols_content[c] #data inside each column
        if type(values[0]) in [int,np.int64,np.int16]:
            data_types[c] = 'K'
        elif type(values[0]) in [np.float16,np.float32,np.float64,np.float]:
            data_types[c] = 'D'
        elif type(values[0]) in [str,np.str_]:
            data_types[c] = '100A'
        elif type(values[0]) is bool:
            data_types[c] = 'L'
        else:
            print("Didn't understand datatype: ",type(values[0]),' for ',
                  c, 'o_O')
            exit()

        # test if it's the ID, in this case convert to integer
        if 'ID' in c:# or 'id' in c: I comment this out because 'confidence' contains id in it
            data_types[c] = '100A'#'K'
            temp = all_cols_content[c] #all MUSE IDs, UV IDs (each kind at a time)
            try:
                all_cols_content[c] = np.asarray(temp,dtype=int)
                sorted_ids = sorted(np.asarray(temp,dtype=int)) #same as temp but sorted in ascending order
                sorted_ids_index = np.argsort(np.asarray(temp,dtype=int))
            except: # in case there are 'XXX_0' (multiple clumps)
                all_cols_content[c] = np.asarray(temp,dtype=str)
                data_types[c] = '20A'
                # sort IDs if str with 'XXX_0'
                ints = []
                for t in temp:
                    try:
                        ints.append(int(t))
                    except:
                        pass
                sorted_ints = sorted(ints)
                sorted_ints_all = []
                for sints in sorted_ints:
                    here_missing = [m for m in temp if str(sints)+'_' in m 
                                    and m[:-2] in str(sints)]
                    sorted_here_missing = sorted(here_missing)
                    sorted_ints_all += [str(sints)]
                    sorted_ints_all += sorted_here_missing
                sorted_ids_index = np.asarray([list(temp).index(ind)
                                               for ind in sorted_ints_all])
                
    cols = []
    cols_sorted = sorted(all_cols_content.keys()) #name of columns sorted in alphabetical order
    for c in cols_sorted:
        if 'ID' in c: # To put ID first
            array = np.asarray(all_cols_content[c])[sorted_ids_index] #same as temp but as array
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)   #column names

    # put RA and DEC next
    for c in cols_sorted:
        if 'RA' in c or 'Ra' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index]
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'DEC' in c or 'Dec' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'Separation' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'lambda' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'z' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'Photometry' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'Confidence' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'Comment' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'ra_noMatch' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
    for c in cols_sorted:
        if 'dec_noMatch' in c:
            array = np.asarray(all_cols_content[c])[sorted_ids_index] 
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
            
    exclude = ['ID','RA','Ra','DEC','Dec','Separation','z','lambda','Photometry','Confidence','Comment','ra_noMatch','dec_noMatch']
    # fill with all other columns, but sorted
    for c in cols_sorted:
        test = [i for i in exclude if i in c]
        if len(test)==0 :
            array = np.asarray(all_cols_content[c])[sorted_ids_index]
            col = fits.Column(name=c,format=data_types[c],array=array)
            cols.append(col)
     
    tbhdu = fits.BinTableHDU.from_columns(cols) #create bin table hdu from scratch
    current_file_name = inspect.getfile(inspect.currentframe()) 
    prihdr = fits.Header()
    prihdr['SCRIPT'] = current_file_name
    prihdr['DATE'] = str(time.strftime("%d/%m/%Y"))
    prihdr['TIME'] = str(time.strftime("%H:%M:%S"))
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu,tbhdu]) #create primary hdu
    thdulist.writeto(f,overwrite=True)
    #hdulist.info() one primary hdu and one bin table hdu
    
    

def modify_output_table(d,f,hdu,position):
    '''
    reads table and modifies it 
    I can pass an existing table to the function and only rewrites specific keys 
    (not the entire table) or appends new stuff if those keys do not exist.
    d = dictionary that is gonna be saved in the output file, the counterparts
    f = name and directory of output file
    position = number of objects we have classified till we called this function
    Check:
    https://thispointer.com/python-how-to-add-append-key-value-pairs-in-dictionary-using-dict-update/#:~:text=We%20can%20add%20%2F%20append%20new,operator%20and%20update()%20function.
    to update, append and write new values to dictionaries
    '''
    
    #read ids+dictionary with key-value pairs (dictionary with IDs as keys, each of which contains another dictionary with the individual values from the columns for each ID)
    info_ids=read_fits_table2dict(f,hdu)[1] 
    #read key-value from counterpart dictionary
    all_cols_content = reshape_dict_in_dict(d,'ID')     
    #I put all_cols_content in the same shape as info_ids
    shaped_all_cols_content = {} 
    for i in range(len(all_cols_content['ID'])):
        id_params = {}
        for c in all_cols_content.keys():
            id_params[c] = all_cols_content[c][i]
        shaped_all_cols_content[str(all_cols_content['ID'][i])] = id_params 
  
    #check if the ID in the counterpart is already in the output table  
    if all_cols_content['ID'][position] in info_ids: #if ID counterpart is already in output file   
        for key, value in all_cols_content.items():  
            #replace the values
            if info_ids[all_cols_content['ID'][position]][key] == all_cols_content[key][position]:
                info_ids[all_cols_content['ID'][position]] = shaped_all_cols_content[all_cols_content['ID'][position]]
            write_table_from_dict_in_dict(info_ids,f)
            print('ID '+str(all_cols_content['ID'][position])+' info replaced in',f)
            break
    else:             
        for i, item in enumerate(shaped_all_cols_content):
            if item not in info_ids.keys():        
                # Append the value in list         
                info_ids.update(shaped_all_cols_content) 
        write_table_from_dict_in_dict(info_ids,f)              
        print('ID info added to',f)    


def invert_dict(d):
    # reverse keys and values
    inv_map = {v: k for k, v in d.items()}
    return inv_map


def reshape_dict_in_dict(old_dict,keyname):
    # make a dictionary of lists from a dictionary of dictionaries
    
    new_vals = list(old_dict) # the old keys
    new_dict = {keyname:new_vals}
    for od in old_dict[new_vals[0]]: # initiate lists
        new_dict[od] = []

    for nv in new_vals:
        for od in old_dict[new_vals[0]]:
            temp = new_dict[od]
            temp.append(old_dict[nv][od])
            new_dict[od] = temp

    return new_dict


def read_image_fits_file(image,hdu):
    
    # open segmentation map
    hdulist = fits.open(image)
    data = hdulist[hdu].data
    header = hdulist[hdu].header
    hdulist.close()

    return data,header  
