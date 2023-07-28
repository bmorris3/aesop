
import numpy as np
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = 1e10   # Otherwise would only show first 50 values

__all__ = ['get_duncan_catalog', 'sindex_catalog', 'query_catalog_for_object',
           'get_k2_epic_catalog']

sindex_catalog = None
k2_epic_table = None
duncan1991 = 'III/159A'
huber2016 = 'J/ApJS/224/2/table5'


def get_duncan_catalog():
    """
    Parameters
    ----------

    Returns
    -------
    """
    global sindex_catalog

    if sindex_catalog is None:
        duncan1991 = 'III/159A'
        sindex_catalog = Vizier(catalog=duncan1991,
                                columns=["*", "Bmag", "Vmag"],
                                row_limit=1e10
                                ).query_constraints()[0]
    return sindex_catalog


def query_catalog_for_object(identifier, catalog=duncan1991):
    """
    Parameters
    ----------
    identifier : str

    catalog : str (optional)

    Returns
    -------

    """
    query = Vizier.query_object(identifier, catalog=catalog)

    if len(query) > 0:
        return query[0][0]
    else:
        return dict(Smean=np.nan, Smin=np.nan, Smax=np.nan)


def get_k2_epic_catalog():
    """
    Parameters
    ----------

    Returns
    -------
    """
    global k2_epic_table

    if k2_epic_table is None:
        catalogs = Vizier.get_catalogs(huber2016)
        k2_epic_table = catalogs[0]  # This is the table with the data
        k2_epic_table.add_index('EPIC')
    return k2_epic_table