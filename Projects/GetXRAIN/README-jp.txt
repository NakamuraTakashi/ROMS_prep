### XバンドMPレーダデータ検索・ダウンロードスクリプト　（2018年11月） ###

    XバンドMPレーダデータの検索・ダウンロードを行うためのPython 3による
  スクリプトです。検索用スクリプトとダウンロード用スクリプトの2つのスク
  リプトにより構成されます。検索用スクリプトによりファイルリストを出力
  し、ダウンロード用スクリプトの引数としてファイル名を指定することでデー
  タのダウンロードをします。

1. データ検索スクリプト - xrain-ls.py
  説明
    指定されたデータ種別、領域・地点のファイル名リストを表示します。
  書式
    xrain-ls.py [-options] {{composite|composite_cs}/AREA|{raw|intermediated}/SITE} ...
  引数
    データ種別名、領域名・地点名を/で区切って指定します。
    データ種別名としては以下のものが指定できます。
      ・composite （合成雨量データ、Xバンドのみ）
      ・raw （rawデータ）
      ・intermediated （一次処理データ）
      ・composite_cx （合成雨量データ、Cバンド+Xバンド）
    領域名・地点名としては以下のものが指定できます。
    （compositeに対して）
      ・SAPPORO
      ・KURIKOMA
      ・FUKUSHIMA
      ・TOUKYOU
      ・NIIGATA
      ・HOKURIKU
      ・SHIZUOKA
      ・NAGOYA
      ・KINKI
      ・OKAYAMA
      ・HIROSHIMA
      ・FUKUOKA
      ・KUMAMOTO
      ・OOSUMI
    （raw/intermediatedに対して）
      ・KITAHIRO
      ・ISHIKARI
      ・ICHIHASAMA
      ・ICHINOSEKI
      ・IWANUMA
      ・WAKUYA
      ・MORIOKA
      ・TAKANOSU
      ・DATE
      ・TAMURA
      ・FUNABASHI
      ・KANTOU
      ・SHINYOKO
      ・UJIIE
      ・YATTAJIMA
      ・KYOUGASE
      ・NAKANOKUTI
      ・MIZUHASHI
      ・NOUMI
      ・FUJINOMIYA
      ・KANUKI
      ・SHIZUKITA
      ・HAMAMATSU
      ・ANJOU
      ・BISAI
      ・SUZUKA
      ・JUUBUSAN
      ・KATSURAGI
      ・ROKKO
      ・TANOKUCHI
      ・KUMAYAMA
      ・TSUNEYAMA
      ・NOGAIBARA
      ・USHIO
      ・FURUTSUKI
      ・KAZASI
      ・KUSENBU
      ・SUGADAKE
      ・UKI
      ・YAMAGA
      ・SAKURAJIMA
    （composite_cxに対して）
      ・HOKAIDO
      ・TOHOKU
      ・KANTO
      ・CYUGOKU
      ・KYUSYU
      ・OKINAWA
  オプション
    -h, --help
        使用法を表示します。
    -f YYYY-MM-DDThh:mm, --from=YYYY-MM-DDThh:mm
        検索対象とする期間の開始時刻（年月日時分）を指定します。※
    -t YYYY-MM-DDThh:mm, --to=YYYY-MM-DDThh:mm
        検索対象とする期間の終了時刻（年月日時分）を指定します。※
    -s MINUTES, --span=MINUTES
        検索対象とする期間の長さを分で指定します。※
    -e ELxxxxxx, --ensemble=ELxxxxxx
        検索対象とする仰角番号を指定します。複数の指定が可能です。未指
        定の場合には全ての仰角番号を検索対象とします。
    -n FILE, --netrc=FILE
        ユーザ設定ファイル名を指定します。ユーザ設定ファイルは.netrcに
        準拠します。machineにはxrain.diasjp.netを指定してください。未指
        定の場合にはユーザのホームディレクトリにある.netrcとなります。
    -u USERNAME, --user=USERNAME
        DIASアカウント名を指定します。
    --version
        スクリプトのバージョンを表示します。

      ※ 検索対象とする期間について
          開始時刻、終了時刻、長さが全て指定された場合、開始時刻、終了
        時刻が優先されます。2つ以上が未指定の場合、終了時刻が未指定の場
        合には現在の時刻、長さが未指定の場合には59分間として検索対象期
        間が設定されます。

2. データダウンロードスクリプト - xrain-dl.py
  説明
    指定されたファイルをダウンロードします。ファイルはtarによりアーカイ
    ブされて出力されます。
  書式
    xrain-dl.py [-options] file ...
  引数
    データ検索スクリプトによって表示されたファイル名を指定します。複数
    のファイル名が指定可能ですが、10000ファイル以上は指定しないでくださ
    い。
  オプション
    -h, --help
        使用法を表示します。
    -o FILE, --output=FILE
        ダウンロードしたファイルを出力するファイル名を指定します。未指
        定の場合には標準出力に出力されます。
    -n FILE, --netrc=FILE
        ユーザ設定ファイル名を指定します。ユーザ設定ファイルは.netrcに
        準拠します。machineにはxrain.diasjp.netを指定してください。未指
        定の場合にはユーザのホームディレクトリにある.netrcとなります。
    -u USERNAME, --user=USERNAME
        DIASアカウント名を指定します。
    --version
        スクリプトのバージョンを表示します。

3. 使用例
  2017年4月1日0:00から2017年4月1日2:59までの北広島（KITAHIRO）および石
  狩（ISHIKARI）のrawデータ（raw）をダウンロードする。
  > xrain-ls.py -f 2017-04-01T00:00 -t 2017-04-01T02:59 raw/{KITAHIRO,ISHIKARI} > list.txt
  > xrain-dl.py -o files.tar `cat list.txt`
