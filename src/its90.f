c ITS90:  1990 temperature scale from the 1968 IPTS68 scale
      function t90(t68)
      t90 = 1.00024*t68
      return
      end
c
c IPTS68:  1968 temperature scale from the 1990 temperature scale
      function t68(t90)
      t68 = 0.99976*t90
      return
      end

