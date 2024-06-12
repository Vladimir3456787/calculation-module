from pymongo import MongoClient
from pymongo.server_api import ServerApi

uri = "mongodb+srv://vovabalaxoncev:Thcvovan7777@cluster0.u499jdc.mongodb.net/?retryWrites=true&w=majority"
class MongoConnection:
    def __init__(self, uri):
        self.client = MongoClient(uri, server_api=ServerApi('1'))
        self.db = self.client['BazaProv']
        self.collection_density = self.db['Rocks']
        self.collection_krepi = self.db['Krepi']
        self.collection_result = self.db['Result']
        self.collection_stratigraphy = self.db['Stratigraphy']
        self.all_stratigraphy = list(self.collection_stratigraphy.find({}, {'_id': 0,'Number': 1, 'X': 1, 'Y': 1}))
        self.collection_areas = self.db['Areas']
        self.collection_analytics = self.db['Analytics']
        self.number_areas = self.collection_areas.distinct('Numbers')
        self.info_areas = self.collection_areas.distinct('Area_info')
        self.collection_areas_full = list(self.collection_areas.find({}, {'_id': 0, 'Area_info': 1, 'Numbers': 1}))
        self.column_X = self.collection_stratigraphy.distinct('X')
        self.column_Y = self.collection_stratigraphy.distinct('Y')
        self.column_number = self.collection_stratigraphy.distinct('Number')
        self.docs = self.collection_density.distinct('Породы')
        self.collection_density_full = list(self.collection_density.find({}, {'_id': 0,'Породы': 1, 'p': 1, 'f': 1, 'p_v': 1 }))
        self.plot = self.collection_density.distinct('Плотность')
        self.name_krep = self.collection_krepi.distinct('Название')
        self.characteristic = self.collection_density.distinct('Характеристики')

    def get_collection_areas(self):
        return self.collection_areas_full

    def get_collection_density(self):
        return self.collection_density_full

    # def get_p(self):
    #     p = [doc['p'] for doc in self.p]
    #     return p
    # def get_f(self):
    #     f = [doc['f'] for doc in self.f]
    #     return f
    def get_collection_analytics(self):
        return self.collection_analytics
    def get_number_areas(self):
        return self.number_areas
    def get_areas(self):
        return self.info_areas
    def get_column_x(self):
        unique_x = [doc['X'] for doc in self.all_stratigraphy]
        return unique_x
    def get_column_y(self):
        unique_y = [doc['Y'] for doc in self.all_stratigraphy]
        return unique_y
    def get_column_number(self):
        return self.column_number
    def get_docs(self):
        return self.docs

    def get_plot(self):
        return self.plot

    def get_name_krep(self):
        return self.name_krep

    def get_density_collection(self):
        return self.collection_density

    def get_krepi_collection(self):
        return self.collection_krepi

    def get_result_collection(self):
        return self.collection_result