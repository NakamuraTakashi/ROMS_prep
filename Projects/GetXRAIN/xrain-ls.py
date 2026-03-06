
import re
import netrc
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import urllib.parse
import http.cookiejar
import html.parser
import optparse
import getpass
import xml.etree.ElementTree

class CASLoginParser(html.parser.HTMLParser):
    def __init__(self):
        html.parser.HTMLParser.__init__(self)
        self.action = None
        self.data = {}

    def handle_starttag(self, tagname, attribute):
        if tagname.lower() == 'form':
            attribute = dict(attribute)
            if 'action' in attribute:
                self.action = attribute['action']
        elif tagname.lower() == 'input':
            attribute = dict(attribute)
            if 'name' in attribute and 'value' in attribute:
                self.data[attribute['name']] = attribute['value']

class DIASAccess():
    def __init__(self, username, password):
        self.__cas_url = 'https://auth.diasjp.net/cas/login?'
        self.__username = username
        self.__password = password
        cj = http.cookiejar.CookieJar()
        self.__opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))

    def open(self, url, data=None):
        response = self.__opener.open(url, data.encode('utf-8'))
        response_url = response.geturl()

        if response_url != url and response_url.startswith(self.__cas_url):
            # redirected to CAS login page
            response = self.__login_cas(response)
            if data != None:
                # If POST (data != None), need reopen
                response.close()
                response = self.__opener.open(url, data.encode('utf-8'))

        return response

    def __login_cas(self, response):
        parser = CASLoginParser()
        parser.feed(response.read().decode('utf-8'))
        parser.close()

        if parser.action == None:
            raise LoginError('Not login page')

        action_url = urllib.parse.urljoin(response.geturl(), parser.action)
        data = parser.data
        data['username'] = self.__username
        data['password'] = self.__password

        response.close()
        response = self.__opener.open(action_url, urllib.parse.urlencode(data).encode('utf-8'))

        if response.geturl() == action_url:
            raise LoginError('Authorization fail')

        return response

class LoginError(Exception):
    def __init__(self, e):
        Exception.__init__(self, e)

if __name__ == '__main__':
    area = {"TOUKYOU":          "TOUKYOU001",
            "HOKURIKU":         "HOKURIKU01",
            "NAGOYA":           "NAGOYA0001",
            "KINKI":            "KINKI00001",
            "KURIKOMA":         "KURIKOMA01",
            "SHIZUOKA":         "SHIZUOKA01",
            "OKAYAMA":          "OKAYAMA001",
            "HIROSHIMA":        "HIROSHIMA1",
            "FUKUOKA":          "FUKUOKA001",
            "OOSUMI":           "OOSUMI0001",
            "NIIGATA":          "NIIGATA001",
            "SAPPORO":          "SAPPORO001",
            "FUKUSHIMA":        "FUKUSHIMA1",
            "KUMAMOTO":         "KUMAMOTO01",
        }
    site = {"KANTOU":           "TOUKYOU001/KANTOU0000",
            "SHINYOKO":         "TOUKYOU001/SHINYOKO00",
            "YATTAJIMA":        "TOUKYOU001/YATTAJIMA0",
            "UJIIE":            "TOUKYOU001/UJIIE00000",
            "FUNABASHI":        "TOUKYOU001/FUNABASHI0",
            "MIZUHASHI":        "HOKURIKU01/MIZUHASHI0",
            "NOUMI":            "HOKURIKU01/NOUMI00000",
            "BISAI":            "NAGOYA0001/BISAI00000",
            "ANJOU":            "NAGOYA0001/ANJOU00000",
            "SUZUKA":           "NAGOYA0001/SUZUKA0000",
            "ROKKO":            "KINKI00001/ROKKO00000",
            "KATSURAGI":        "KINKI00001/KATSURAGI0",
            "JUUBUSAN":         "KINKI00001/JUUBUSAN00",
            "TANOKUCHI":        "KINKI00001/TANOKUCHI0",
            "WAKUYA":           "KURIKOMA01/WAKUYA0000",
            "IWANUMA":          "KURIKOMA01/IWANUMA000",
            "MORIOKA":          "KURIKOMA01/MORIOKA000",
            "TAKANOSU":         "KURIKOMA01/TAKANOSU00",
            "ICHINOSEKI":       "KURIKOMA01/ICHINOSEKI",
            "ICHIHASAMA":       "KURIKOMA01/ICHIHASAMA",
            "SHIZUKITA":        "SHIZUOKA01/SHIZUKITA0",
            "KANUKI":           "SHIZUOKA01/KANUKI0000",
            "FUJINOMIYA":       "SHIZUOKA01/FUJINOMIYA",
            "HAMAMATSU":        "SHIZUOKA01/HAMAMATSU0",
            "KUMAYAMA":         "OKAYAMA001/KUMAYAMA00",
            "TSUNEYAMA":        "OKAYAMA001/TSUNEYAMA0",
            "NOGAIBARA":        "HIROSHIMA1/NOGAIBARA0",
            "USHIO":            "HIROSHIMA1/USHIO00000",
            "KUSENBU":          "FUKUOKA001/KUSENBU000",
            "SUGADAKE":         "FUKUOKA001/SUGADAKE00",
            "FURUTSUKI":        "FUKUOKA001/FURUTSUKI0",
            "KAZASI":           "FUKUOKA001/KAZASI0000",
            "SAKURAJIMA":       "OOSUMI0001/SAKURAJIMA",
            "KYOUGASE":         "NIIGATA001/KYOUGASE00",
            "NAKANOKUTI":       "NIIGATA001/NAKANOKUTI",
            "KITAHIRO":         "SAPPORO001/KITAHIRO00",
            "ISHIKARI":         "SAPPORO001/ISHIKARI00",
            "DATE":             "FUKUSHIMA1/DATE000000",
            "TAMURA":           "FUKUSHIMA1/TAMURA0000",
            "YAMAGA":           "KUMAMOTO01/YAMAGA0000",
            "UKI":              "KUMAMOTO01/UKI0000000",
            }
    area_cx = {"HOKAIDO":       "HOKAIDO001",
               "TOHOKU":        "TOHOKU0001",
               "KANTO":         "KANTO00001",
               "CYUGOKU":       "CYUGOKU001",
               "KYUSYU":        "KYUSYU0001",
               "OKINAWA":       "OKINAWA001",
               }
    target = {'composite':      {'root':        'composite/',
                                 'location':    area},
              'raw':            {'root':        'kyoku1/',
                                 'location':    site},
              'intermediated':  {'root':        'kyoku2/',
                                 'location':    site},
              'composite_cx':   {'root':        'composite_cx/',
                                 'location':    area_cx}}
    host = 'xrain.diasjp.net'
    url = 'http://' + host + '/original/search.cgi'

    usage ='''usage: %prog [options] target ...
  target\t{composite|composite_cx}/AREA or {raw|intermediated}/SITE'''
    version='%prog 18.1116'
    parser = optparse.OptionParser(usage=usage, version=version)
    parser.add_option('-f', '--from', dest='start',
                      help='specify the start time of observation period',
                      metavar='YYYY-MM-DDThh:mm')
    parser.add_option('-t', '--to', dest='end',
                      help='specify the end time of observation period', 
                      metavar='YYYY-MM-DDThh:mm')
    parser.add_option('-s', '--span',
                      help='specify the span of observation period',
                      metavar='MINUTES')
    parser.add_option('-e', '--elevation', action='append',
                      help='specify the elevation (EL000000-EL120000)', 
                      metavar='ELxxxxxx')
    parser.add_option('-n', '--netrc', default=None,
                      help='specify the netrc file', metavar='FILE')
    parser.add_option('-u', '--user', default=None,
                      help='specify the DIAS account name',
                      metavar='USERNAME')

    (options, args) = parser.parse_args()

    if len(args) <= 0:
        parser.error('No target specified')

    if options.elevation is None:
        options.elevation = ['EL%02d0000' % i for i in range(13)]

    param = []
    if options.start is not None:
        if not re.match(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}', options.start):
            parser.error('Illegal time format: ' + options.start)
        param.append('start=' + options.start + ':00')
    if options.end is not None:
        if not re.match(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}', options.end):
            parser.error('Illegal time format: ' + options.end)
        param.append('end=' + options.end + ':00')
    if options.span is not None:
        param.append('span=' + options.span)
    for el in options.elevation:
        param.append('elevation[]=' + el)

    for dt in args:
        (tp, lc) = dt.split('/')
        if tp not in target:
            print(tp)
            parser.error('Unknown target: ' + dt)
        if lc not in target[tp]['location']:
            print(lc)
            parser.error('Unknown target: ' + dt)
        param.append('data[]=' + target[tp]['root'] + 
                     target[tp]['location'][lc])

    (login, password) = (None, None)

    try:
        auth = netrc.netrc(options.netrc).authenticators(host)
        if auth is not None:
            (login, account, password) = auth
    except (IOError):
        pass

    if options.user is not None:
        login = options.user
        password = None

    if login is None:
        login = input('Username: ')

    if password is None:
        password = getpass.getpass('Password: ')

    access = DIASAccess(login, password)
    response = access.open(url, '&'.join(param))
    root = xml.etree.ElementTree.fromstring(response.read())
    response.close()

    lst = root.findall('item')
    for i in lst:
        print(i.attrib['id'])
