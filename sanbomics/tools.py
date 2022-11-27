import pandas as pd
import pkg_resources



class id_map(object):
    '''
    initialize the id mapping object
    
    species : string
        required. Name of species: human or mouse
        
    key: string (default: ensembl)
        The id type you are using as input keys.
        Options: symbol, ensembl, entrez
    
    target : string (default: symbol)
        The id type you want returned.
        Options: symbol, ensembl, entrez
        
    '''
    
    
    def __init__(self, species = None, key = 'ensembl', target = 'symbol'):
        
        self.key = key
        
        self.target = target
        
        self.species = species
        
        if self.species not in ['human', 'mouse']:
            raise Exception("set species to one of the following: human, mouse")
        
        
        if self.species == 'human':
        	stream = pkg_resources.resource_stream(__name__, 'data/human_ids.csv')
            self.dataframe = pd.read_csv(stream)
            self.dataframe = self.dataframe.rename(columns={'hgnc_symbol': 'symbol',
                                                            'ensembl_gene_id':'ensembl',
                                                            'entrezgene_id':'entrez'})
        else:
        	stream = pkg_resources.resource_stream(__name__, 'data/mouse_ids.csv')
            self.dataframe = pd.read_csv(stream)
            self.dataframe = self.dataframe.rename(columns={'uniprot_gn_symbol': 'symbol',
                                                            'ensembl_gene_id':'ensembl',
                                                            'entrezgene_id':'entrez'})
        
        
        self.mapper = self.dataframe[(self.dataframe[self.key].notna()) &\
                                     (self.dataframe[self.target].notna())]
        
        
        if target != 'entrez':
            self.mapper = dict(zip(self.mapper[self.key], self.mapper[self.target]))
        else: #because entrez is imported as floats
            self.mapper = dict(zip(self.mapper[self.key],
                                   self.mapper[self.target].map(lambda x: str(round(x)))))
            
            
    def map_column(self, df, column = None):
        '''
        map your mapping object to a pandas dataframe column.
        
        Returns a dataframe with an appended column (target + _map) containing mapped ids.
        
        df: pandas.DataFrame
        
        column: string
            Name of column containing the keys
        
        '''
        df[self.target + '_map'] = df[column].map(self.mapper)
        return df
    
    def map_list(self, l):
        '''
        Returns a mapped list with 'NA' for unmappable keys
        
        l: list
        
        '''
        def mini_mapper(x):
            try:
                return self.mapper[x]
            except:
                return 'NA'
        return [mini_mapper(x) for x in l]





