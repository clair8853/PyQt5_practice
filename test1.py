from PIL import Image, ImageOps
# 2色の画像を混ぜる

# 画像ファイルのパス
image_path1 = "C:/Users/visiu/Documents/Figure_1.png"  # 画像1のパス（Redチャネルとして使用）
image_path2 = 'C:/Users/visiu/Documents/Figure_2.png'  # 画像2のパス（Greenチャネルとして使用）

# 画像を読み込む
img1 = Image.open(image_path1)
img2 = Image.open(image_path2)

# 画像を白黒に変換
bw_img1 = img1.convert('L')
bw_img2 = img2.convert('L')

# 画像を反転させる
inverted_img1 = ImageOps.invert(bw_img1)
inverted_img2 = ImageOps.invert(bw_img2)

# 反転させた白黒画像をそのままRGBに変換
# この時、RGBに変換してチャンネルを扱いやすくする
red_channel = inverted_img1.convert('RGB').split()[0]  # 反転させた画像をRedチャネルに
green_channel = inverted_img2.convert('RGB').split()[1]  # 反転させた画像をGreenチャネルに

# 空のBlueチャネルを作成（全てのピクセルが0の画像）
blue_channel = Image.new('L', bw_img1.size, 0)

# 新しい画像を作成する
merged_image = Image.merge('RGB', (red_channel, green_channel, blue_channel))

# 解像度を10倍に拡大
new_size = (int(merged_image.width * 10), int(merged_image.height * 10))  
high_res_image = merged_image.resize(new_size, Image.Resampling.LANCZOS)

# 高解像度で保存
#high_res_image.save('merged_image.jpg', quality=100)

#Redチャネルの画像保存
red_channel = inverted_img1.convert('RGB').split()[0]  # 反転した画像のRedチャネル
green_channel = Image.new('L', inverted_img1.size, 0)  # 空のGreenチャネル
blue_channel = Image.new('L', inverted_img1.size, 0)   # 空のBlueチャネル

# 新しい画像を作成する
red_image = Image.merge('RGB', (red_channel, green_channel, blue_channel))

# 解像度を10倍に拡大
new_size = (int(red_image.width * 10), int(red_image.height * 10))  
high_res_image1 = red_image.resize(new_size, Image.Resampling.LANCZOS)

# 高解像度で保存



high_res_image.show()


high_res_image1.show()

