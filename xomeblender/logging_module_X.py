import logging
import traceback
from xomeblender.colored_log import ColoredFormatter

# Create top level logger
log = logging.getLogger("InXalizer")

# Add console handler using our custom ColoredFormatter
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
cf = ColoredFormatter("[%(name)s][%(levelname)s] - %(message)s") #(%(filename)s:%(lineno)d)
ch.setFormatter(cf)
log.addHandler(ch)

# Add file handler
fh = logging.FileHandler('InXalizer.log')
fh.setLevel(logging.DEBUG)
ff = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(ff)
log.addHandler(fh)

# Set log level
log.setLevel(logging.INFO)
