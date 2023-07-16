## 使用方法

0. 点击 `Code` - `Download ZIP` 下载压缩包，解压

### 包

1. 将文件移动至 `%APPDATA%/Mathematica/Autoload/Mirion`
2. 重启 Mathematica

- `Kernel/init.m` - 作用是导入 `main.m`
- `main.m` - 本体，里面只有格式化后的代码
- `dev.nb` - 开发文件，附有注释

### 样式美化

1. 将 `DefaultStyleDefinitions` 项修改为 `style.nb` 的位置
2. 执行 `<< Mirion/tf.m`

- `style.nb` - 样式表
- `tf.m` - 传统格式设置脚本
