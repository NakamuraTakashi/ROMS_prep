
import netrc
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import urllib.parse
import http.cookiejar
import html.parser
import optparse
import getpass
import sys

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
    host = 'xrain.diasjp.net'
    url = 'http://' + host + '/original/download.cgi'

    usage ='''usage: %prog [options] file ...'''
    version='%prog 18.0816'
    parser = optparse.OptionParser(usage=usage, version=version)
    parser.add_option('-o', '--output',
                      help='specify the output file', 
                      metavar='FILE')
    parser.add_option('-n', '--netrc', default=None,
                      help='specify the netrc file', metavar='FILE')
    parser.add_option('-u', '--user', default=None,
                      help='specify the DIAS account name',
                      metavar='USERNAME')

    (options, args) = parser.parse_args()

    if len(args) <= 0:
        parser.error('No file specified')
    if len(args) > 10000:
        parser.error('Too many files')

    param = ['archiver=tar']

    for dt in args:
        param.append('file=' + dt)

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

    if options.output is not None:
        f = open(options.output, 'wb')
    else:
        f = sys.stdout

    while True:
        buf = response.read(32768)
        if not buf:
            break

        f.write(buf)

    if options.output is not None:
        f.close()

    response.close()
