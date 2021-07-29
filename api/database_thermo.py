from api.config import Config
from api.thermodynamics import MINE_thermo

mine_thermo = MINE_thermo(mongo_uri=Config.MONGO_URI,
                          postgres_uri=Config.POSTGRES_URI)
