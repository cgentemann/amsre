"""
provides various function for sea surface temperature processing and analysis
"""
import numpy


COASTAL_FLAG_BIT = 6
ICE_FLAG_BIT = 4
CLOUDY_FLAG_BIT = 2
RAINY_FLAG_BIT = 3


def get_ghrsst_cloud_mask(feature):
    '''
    Return the cloud mask of a GHRSST product

    Args:
        feature(cerbere object): Swath or Grid instanciated with a GHRSST dataset from which to extract the cloud mask
    '''
#    if 'NAVO-L2P' in feature._dataStorage.get_collection_id():
        # in NAVOCEANO products, rejection_flag is not filled. We have to use the sst quality level instead.
#        qlt = feature.get_values('proximity_confidence')
#         cloudy = numpy.ma.masked_equal(qlt,2).mask
#        print qlt.min(), qlt.max()
#    else:
    qlt = feature.get_values('rejection_flag')
    print qlt.min(), qlt.max()
    temp = numpy.ma.logical_and( (qlt & (1<<COASTAL_FLAG_BIT)==0), numpy.ma.logical_or((qlt & (1<<CLOUDY_FLAG_BIT)>0),(qlt & (1<<RAINY_FLAG_BIT)>0)) )
    cloudy = numpy.ma.array(temp, mask=numpy.ma.make_mask(temp), copy=0)
    return cloudy

def get_ghrsst_rain_mask(feature):
    '''
    Return the rain mask of a GHRSST product

    Args:
        feature (cerbere object): Swath or Grid instanciated with a GHRSST dataset from which to extract the rain mask
    '''
    qlt = feature.get_values('rejection_flag')
    temp = numpy.ma.logical_and( (qlt & (1<<COASTAL_FLAG_BIT)==0), numpy.ma.logical_or((qlt & (1<<CLOUDY_FLAG_BIT)>0),(qlt & (1<<RAINY_FLAG_BIT)>0)) )
    cloudy = numpy.ma.array(temp, mask=numpy.ma.make_mask(temp), copy=0)
    return cloudy


def get_ghrsst_land_mask(feature):
    '''
    Return the land mask of a GHRSST product

    Args:
        feature(cerbere object): Swath or Grid instanciated with a GHRSST dataset from which to extract the land mask
    '''
    if feature.has_field('rejection_flag'):
        qlt = feature.get_values('rejection_flag')
        temp = (qlt & (1<<COASTAL_FLAG_BIT) > 0)
    else:
        qlt = feature.get_values('l2p_flags')
        LAND_FLAG = 1  #2s
        temp = (qlt & (1<<LAND_FLAG) > 0)
        
    land = numpy.ma.array(temp, mask=numpy.ma.make_mask(temp), copy=0)
    
    return land


def get_ghrsst_ice_mask(feature):
    '''
    Return the ice mask of a GHRSST product

    Args:
        feature(cerbere object): Swath or Grid instanciated with a GHRSST dataset from which to extract the ice mask
    '''
    if feature.has_field('rejection_flag'):
        qlt = feature.get_values('rejection_flag')
        temp = (qlt & (1<<ICE_FLAG_BIT) > 0)
    else:
        qlt = feature.get_values('l2p_flags')
        ICE_FLAG = 2  #4s
        temp = (qlt & (1<<ICE_FLAG) > 0)
    
    ice = numpy.ma.array(temp, mask=numpy.ma.make_mask(temp), copy=0)
    
    return ice


def get_ghrsst_quality_sst( sst, quality, quality_min_level=3, quality_max_level=None):
    '''
    Return sst field filtered from low quality measurements

    Args:
        sst (float): masked array of floats GHRSST dataset from which to extract quality data
        quality (masked_array of floats) :quality level values (proximity_confidence in GDS v1, quality_level in GDS v2)
        quality_min_level (integer): minimum quality level of the data to select
        quality_min_level (integer): maximum quality level of the data to select
    '''
    if quality_max_level is None:
        return numpy.ma.masked_where(quality < quality_min_level, sst)
    elif quality_min_level == quality_max_level:
        return numpy.ma.masked_where(quality != quality_min_level, sst)
    else:
        return numpy.ma.masked_where((quality < quality_min_level) & (quality > quality_max_level), sst)

def get_percentiles(data, percentiles):
    '''
    Return percentiles of a data distribution

    Args:
        data (masked array of floats):data from which to compute the percentiles
        percentiles (array of integer):values of percentile to compute
    '''
    tmp = numpy.sort(data.compressed())
    res = []
    if (len(tmp) > 0):
        for p in percentiles:
            res.append(tmp[int(len(tmp) * (p / 100.0))])
    return res
